#### Molecule Align ####
# Align molecules using the Kabsch algorithm.
# https://en.wikipedia.org/wiki/Kabsch_algorithm

from ovito.data import DataCollection, DataTable
from ovito.modifiers import AffineTransformationModifier
import numpy as np
from ovito.pipeline import ModifierInterface
from traits.api import Bool, Int


class MoleculeAlign(ModifierInterface):
    only_selected = Bool(True, label="Use only selected particles")
    reference_frame = Int(0, label="Reference frame")

    def get_suffix(self, data):
        max_count = -1
        for key in data.attributes:
            if key.startswith("MoleculeAlign"):
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

    def input_caching_hints(self, frame, **kwargs):
        return [self.reference_frame, frame]

    def modify(
        self,
        data: DataCollection,
        frame: int,
        input_slots: dict[str, ModifierInterface.InputSlot],
        data_cache: DataCollection,
        **kwargs,
    ):
        # Add identifier if not present
        if "Particle Identifier" not in data.particles:
            data.particles_.create_property(
                "Particle Identifier", data=np.arange(data.particles.count)
            )
        data_ref = input_slots["upstream"].compute(self.reference_frame)
        if "Particle Identifier" not in data_ref.particles:
            data_ref.particles_.create_property(
                "Particle Identifier", data=np.arange(data_ref.particles.count)
            )

        # Get selections
        if self.only_selected and "Selection" not in data.particles:
            raise ValueError("No selection available. Please define a selection.")

        if self.only_selected:
            selection = data.particles["Selection"] == 1
        else:
            selection = ...

        if self.only_selected:
            selection_ref = data_ref.particles["Selection"] == 1
        else:
            selection_ref = ...

        # get reference points
        pos_ref = data_ref.particles["Position"][selection_ref]
        idx_ref = data_ref.particles["Particle Identifier"][selection_ref]
        pos_ref = pos_ref[np.argsort(idx_ref)]

        # get current points
        pos = data.particles["Position"][selection]
        idx = data.particles["Particle Identifier"][selection]
        pos = pos[np.argsort(idx)]

        # remove translation
        pos_ref_bar = np.mean(pos_ref, axis=0)
        pos_ref -= pos_ref_bar
        pos_bar = np.mean(pos, axis=0)
        pos -= pos_bar

        # Compute covariance matrix
        H = np.dot(pos.T, pos_ref)
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
        pos_ref = data_ref.particles["Position"][selection_ref]
        pos_ref = pos_ref[np.argsort(idx_ref)]

        pos = data.particles["Position"][selection]
        pos = pos[np.argsort(idx)]

        prop_name = f"MoleculeAlign{self.get_suffix(data)}"

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

        # Save RMSD
        if f"{prop_name}.RMSD.array" not in data_cache.attributes:
            data_cache.attributes[f"{prop_name}.RMSD.array"] = np.empty(
                input_slots["upstream"].num_frames
            )
            data_cache.attributes[f"{prop_name}.RMSD.array"][:] = np.nan

            data_cache.attributes[f"{prop_name}.RMSD_prev.array"] = np.empty(
                input_slots["upstream"].num_frames
            )
            data_cache.attributes[f"{prop_name}.RMSD_prev.array"][:] = np.nan

            data_cache.attributes[f"{prop_name}.RMSD_all.array"] = np.empty(
                input_slots["upstream"].num_frames
            )
            data_cache.attributes[f"{prop_name}.RMSD_all.array"][:] = np.nan

        rmsd_array = data_cache.attributes[f"{prop_name}.RMSD.array"]
        rmsd_array[frame] = rmsd
        table = data.tables.create(
            identifier=f"{prop_name}.RMSD",
            plot_mode=DataTable.PlotMode.Scatter,
            title="MoleculeAlign RMSD",
        )
        table.x = table.create_property(
            "Frame", data=np.arange(input_slots["upstream"].num_frames)
        )
        table.y = table.create_property("RMSD", data=rmsd_array)

        rmsd_array_all = data_cache.attributes[f"{prop_name}.RMSD_all.array"]
        rmsd_array_all[frame] = rmsd_all
        table = data.tables.create(
            identifier=f"{prop_name}.RMSD_all",
            plot_mode=DataTable.PlotMode.Scatter,
            title="MoleculeAlign RMSD All",
        )
        table.x = table.create_property(
            "Frame", data=np.arange(input_slots["upstream"].num_frames)
        )
        table.y = table.create_property("RMSD All", data=rmsd_array_all)
