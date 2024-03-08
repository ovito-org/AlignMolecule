#### Align Molecule ####
# Align a molecule using Kabsch algorithm.
# https://en.wikipedia.org/wiki/Kabsch_algorithm

import platform

if tuple(int(v) for v in platform.python_version_tuple()) <= (3, 8):
    from __future__ import annotations  # noqa: F404

from ovito.data import DataCollection
from ovito.modifiers import AffineTransformationModifier
import numpy as np
from ovito.pipeline import ModifierInterface
from traits.api import Bool, Int


class AlignMolecule(ModifierInterface):
    only_selected = Bool(True, label="Use only selected particles")
    reference_frame = Int(0, label="Reference frame")

    def get_suffix(self, data):
        max_count = -1
        for key in data.attributes:
            if key.startswith("AlignMolecule"):
                key = key.split(".")
                if len(key) == 1:
                    max_count = 0
                else:
                    count = int(key[-1])
                    if count > max_count:
                        max_count = count
        if max_count != -1:
            return f".{max_count+1}"
        return ""

    def get_selection(self, col: DataCollection):
        if "Particle Identifier" in col.particles:
            if self.only_selected:
                selection = np.where(col.particles["Selection"] != 0)[0]
                return selection[
                    np.argsort(col.particles["Particle Identifier"][selection])
                ]
            else:
                return np.argsort(col.particles["Particle Identifier"])
        else:
            if self.only_selected:
                return np.where(col.particles["Selection"] != 0)[0]
            else:
                return ...

    def input_caching_hints(self, frame, **kwargs):
        return [self.reference_frame, frame]

    def modify(
        self,
        data: DataCollection,
        input_slots: dict[str, ModifierInterface.InputSlot],
        **kwargs,
    ):
        # Get selections
        if self.only_selected and "Selection" not in data.particles:
            raise ValueError("No selection available. Please define a selection.")

        selection = self.get_selection(data)
        data_ref = input_slots["upstream"].compute(self.reference_frame)
        selection_ref = self.get_selection(data_ref)

        # get reference points
        pos_ref = data_ref.particles["Position"][selection_ref]

        # get current points
        pos = data.particles["Position"][selection]

        # calculate translation
        pos_ref_bar = np.mean(pos_ref, axis=0)
        pos_bar = np.mean(pos, axis=0)

        # Compute covariance matrix
        H = np.dot((pos - pos_bar).T, (pos_ref - pos_ref_bar))
        U, S, Vt = np.linalg.svd(H)

        # Compute rotation matrix
        d = np.sign(np.linalg.det(Vt) * np.linalg.det(U))
        R = np.dot(
            np.dot(
                Vt.T,
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]]),
            ),
            U.T,
        )

        # Apply rotation to current points
        transform = np.zeros((3, 4))
        transform[:3, :3] = R
        data.apply(AffineTransformationModifier(transformation=transform))

        # Translate points to reference position
        pos = data.particles["Position"][selection]
        translate = pos_ref_bar - np.mean(pos, axis=0)

        transform = np.zeros((3, 4))
        np.fill_diagonal(transform, 1)
        transform[:, 3] = translate
        data.apply(AffineTransformationModifier(transformation=transform))

        # RMSD selection
        pos = data.particles["Position"][selection]
        pos_ref = data_ref.particles["Position"][selection_ref]

        prop_name = f"AlignMolecule{self.get_suffix(data)}"

        rmsd = np.mean(np.square(pos_ref - pos))
        data.attributes[f"{prop_name}.RMSD"] = rmsd

        # RMSD all
        mapping = data_ref.particles.remap_indices(data.particles)
        pos = data.particles["Position"]
        pos_ref = data_ref.particles["Position"][mapping]

        rmsd_all = np.mean(np.square(pos_ref - pos), axis=1)
        data.particles_.create_property("RMSD", data=rmsd_all)
        rmsd_all = np.mean(rmsd_all)
        data.attributes[f"{prop_name}.RMSD_all"] = rmsd_all
