import numpy as np
import pytest
from ovito.io import import_file
from ovito.modifiers import ExpressionSelectionModifier, UnwrapTrajectoriesModifier

from AlignMolecule import AlignMolecule


@pytest.fixture()
def setup_pipeline():
    pipeline = import_file("RDX.reax.dump")
    pipeline.modifiers.append(UnwrapTrajectoriesModifier())

    return pipeline, "||".join(
        (f"ParticleIdentifier == {i}" for i in (5, 4, 6, 1, 2, 3))
    )


def test_rmsd(setup_pipeline):
    pipeline, expression = setup_pipeline
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=expression))
    pipeline.modifiers.append(AlignMolecule())
    data = pipeline.compute(6)
    assert np.isclose(data.attributes["AlignMolecule.RMSD"], 0.004767366293272777)
    assert np.isclose(data.attributes["AlignMolecule.RMSD_all"], 16.719164645542197)


def test_rmsd_per_particle(setup_pipeline):
    pipeline, expression = setup_pipeline
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=expression))
    pipeline.modifiers.append(AlignMolecule())
    data = pipeline.compute(6)
    ref_rmsd = np.array(
        (
            0.000118518,
            0.00419554,
            0.0135968,
            0.00891107,
            0.000861882,
            0.00092043,
            0.0394605,
            0.130131,
            0.000545504,
            1.62779,
            1.65908,
            0.735138,
            0.0742365,
            1.31378,
            1.68638,
            0.0229692,
            0.00875083,
            0.141561,
            0.0884447,
            0.00665578,
            0.0328006,
            49.7822,
            50.1849,
            47.8836,
            45.5941,
            14.2866,
            16.1272,
            54.9405,
            11.949,
            17.8095,
            57.5626,
            7.5454,
            11.1038,
            14.1957,
            15.3778,
            23.7397,
            50.5864,
            48.5399,
            46.6396,
            50.2729,
            45.7907,
            15.6758,
            16.1709,
            27.2426,
            33.0774,
            25.4604,
            16.4809,
            12.7509,
            32.4243,
            25.1584,
            9.6679,
            28.816,
            42.7217,
            36.3903,
            16.2667,
            3.89775,
            16.9398,
            15.5234,
            13.616,
            38.7949,
            42.375,
            19.1103,
            12.7872,
            8.95033,
            6.5917,
            10.7216,
            20.3724,
            22.553,
            17.9234,
            5.32852,
            28.1354,
            21.3825,
            5.38776,
            6.00336,
            26.3668,
            41.7576,
            25.1697,
            22.5754,
            5.87788,
            10.1258,
            10.2545,
            10.3131,
            26.2214,
            25.3171,
            6.76095,
            7.26998,
            8.62313,
            10.7311,
            9.59805,
            7.50536,
            8.62286,
            14.3764,
            9.02929,
            16.1702,
            5.09519,
            20.16,
            13.7649,
            12.7437,
            7.48286,
            6.73011,
            5.8531,
            8.47722,
            9.64693,
            9.31856,
            11.3662,
        )
    )

    for i in range(data.particles.count):
        ref_val = ref_rmsd[i]
        idx = np.where(data.particles["Particle Identifier"] == i + 1)[0]
        assert len(idx) == 1
        assert np.isclose(data.particles["RMSD"][idx], ref_val)


def test_rmsd_all(setup_pipeline):
    pipeline, _ = setup_pipeline
    pipeline.modifiers.append(AlignMolecule(only_selected=False))
    data = pipeline.compute(6)
    assert np.isclose(data.attributes["AlignMolecule.RMSD"], 2.11002521330366)
    assert np.isclose(data.attributes["AlignMolecule.RMSD_all"], 2.11002521330366)
