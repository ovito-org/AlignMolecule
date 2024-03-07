from ovito.data import DataCollection, DataTable
from ovito.modifiers import AffineTransformationModifier
import numpy as np
from ovito.pipeline import ModifierInterface


class MoleculeAlign(ModifierInterface):
    # https://en.wikipedia.org/wiki/Kabsch_algorithm

    def input_caching_hints(self, frame, input_slots, **kwargs):
        return [0, frame]

    def modify(
        self,
        data: DataCollection,
        frame: int,
        input_slots: dict[str, ModifierInterface.InputSlot],
        data_cache: DataCollection,
        **kwargs,
    ):
        assert "Selection" in data.particles

        # get reference points
        data_ref = input_slots["upstream"].compute(0)
        pos_ref = data_ref.particles["Position"][data_ref.particles["Selection"] == 1]
        idx_ref = data_ref.particles["Particle Identifier"][
            data_ref.particles["Selection"] == 1
        ]
        pos_ref = pos_ref[np.argsort(idx_ref)]

        # get current points
        pos = data.particles["Position"][data.particles["Selection"] == 1]
        idx = data.particles["Particle Identifier"][data.particles["Selection"] == 1]
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
        pos = data.particles["Position"][data.particles["Selection"] == 1]
        translate = pos_ref_bar - np.mean(pos, axis=0)

        transform = np.zeros((3, 4))
        np.fill_diagonal(transform, 1)
        transform[:, 3] = translate
        data.apply(AffineTransformationModifier(transformation=transform))

        # RMSD selection
        pos_ref = data_ref.particles["Position"][data_ref.particles["Selection"] == 1]
        pos_ref = pos_ref[np.argsort(idx_ref)]

        pos = data.particles["Position"][data.particles["Selection"] == 1]
        pos = pos[np.argsort(idx)]

        rmsd = np.mean(np.square(pos_ref - pos))
        data.attributes["PointAlignment.RMSD"] = rmsd

        # RMSD all
        pos = data.particles["Position"][
            np.argsort(data.particles["Particle Identifier"])
        ]
        pos_ref = data_ref.particles["Position"][
            np.argsort(data_ref.particles["Particle Identifier"])
        ]

        rmsd_all = np.mean(np.square(pos_ref - pos))
        data.attributes["PointAlignment.RMSD_all"] = rmsd_all

        # Save RMSD
        if "PointAlignment.RMSD.array" not in data_cache.attributes:
            data_cache.attributes["PointAlignment.RMSD.array"] = np.empty(
                input_slots["upstream"].num_frames
            )
            data_cache.attributes["PointAlignment.RMSD.array"][:] = np.nan

            data_cache.attributes["PointAlignment.RMSD_prev.array"] = np.empty(
                input_slots["upstream"].num_frames
            )
            data_cache.attributes["PointAlignment.RMSD_prev.array"][:] = np.nan

            data_cache.attributes["PointAlignment.RMSD_all.array"] = np.empty(
                input_slots["upstream"].num_frames
            )
            data_cache.attributes["PointAlignment.RMSD_all.array"][:] = np.nan

        rmsd_array = data_cache.attributes["PointAlignment.RMSD.array"]
        rmsd_array[frame] = rmsd
        table = data.tables.create(
            identifier="PointAlignment.RMSD",
            plot_mode=DataTable.PlotMode.Line,
            title="PointAlignment RMSD",
        )
        table.x = table.create_property(
            "Frame", data=np.arange(input_slots["upstream"].num_frames)
        )
        table.y = table.create_property("RMSD", data=rmsd_array)

        rmsd_array_all = data_cache.attributes["PointAlignment.RMSD_all.array"]
        rmsd_array_all[frame] = rmsd_all
        table = data.tables.create(
            identifier="PointAlignment.RMSD_all",
            plot_mode=DataTable.PlotMode.Line,
            title="PointAlignment RMSD All",
        )
        table.x = table.create_property(
            "Frame", data=np.arange(input_slots["upstream"].num_frames)
        )
        table.y = table.create_property("RMSD All", data=rmsd_array_all)
