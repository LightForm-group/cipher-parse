from textwrap import dedent
import pytest

import numpy as np
from cipher_parse.cipher_input import (
    CIPHERGeometry,
    CIPHERInput,
    InterfaceDefinition,
    MaterialDefinition,
    PhaseTypeDefinition,
    compress_1D_array_string,
    decompress_1D_array_string,
)
from cipher_parse.discrete_voronoi import DiscreteVoronoi
from cipher_parse.errors import (
    GeometryDuplicateMaterialNameError,
    GeometryExcessTargetVolumeFractionError,
    GeometryMissingPhaseAssignmentError,
    GeometryNonUnitTargetVolumeFractionError,
    GeometryVoxelPhaseError,
    MaterialPhaseTypeFractionError,
    MaterialPhaseTypeLabelError,
    MaterialPhaseTypePhasesMissingError,
)


def test_expected_compresss_1D_array():
    """Test compression as expected."""
    arr = np.array([1, 1, 1, 1, 2, 2, 1, 2, 3, 1, 3, 3, 2, 2, 2, 1, 1, 4])
    arr_str = compress_1D_array_string(arr)
    expected = dedent(
        """
        4 of 1
        2 of 2
        1
        2
        3
        1
        2 of 3
        3 of 2
        2 of 1
        4
        """
    ).strip()
    assert arr_str == expected


def test_round_trip_1D_array():
    """Test round-trip compress/decompress of 1D array."""
    arr = np.random.choice(3, size=100)  # likely to get consecutive repeats
    arr_str = compress_1D_array_string(arr)
    arr_reload = decompress_1D_array_string(arr_str)
    assert np.all(arr_reload == arr)


def get_boiler_plate_geometry_args(
    size=[1, 1], grid_size=[128, 128], num_phases=10, interfaces=None
):
    voronoi_obj = DiscreteVoronoi.from_random(
        size=size, grid_size=grid_size, num_regions=num_phases
    )
    out = {
        "voxel_phase": voronoi_obj.region_ID,
        "size": size,
        "interfaces": interfaces or [],
    }
    return out


def test_material_definition_raise_on_phase_type_same_type_label():
    with pytest.raises(MaterialPhaseTypeLabelError):
        MaterialDefinition(
            name="mat1",
            properties={},
            phase_types=[
                PhaseTypeDefinition(type_label="1"),
                PhaseTypeDefinition(type_label="1"),
            ],
        )


def test_material_definition_raise_on_no_remaining_type_fraction():
    with pytest.raises(MaterialPhaseTypeFractionError):
        MaterialDefinition(
            name="mat1",
            properties={},
            phase_types=[
                PhaseTypeDefinition(type_label="1", target_type_fraction=1.0),
                PhaseTypeDefinition(type_label="2"),
            ],
        )


def test_material_definition_raise_on_type_fraction_over_one():
    with pytest.raises(MaterialPhaseTypeFractionError):
        MaterialDefinition(
            name="mat1",
            properties={},
            phase_types=[
                PhaseTypeDefinition(type_label="1", target_type_fraction=0.3),
                PhaseTypeDefinition(type_label="2", target_type_fraction=1.0),
            ],
        )


def test_geometry_raise_on_missing_phases_second_material():
    materials = [
        MaterialDefinition(
            name="mat1",
            properties={},
            phases=[1, 2, 3],
        ),
        MaterialDefinition(
            name="mat2",
            properties={},
        ),
    ]
    with pytest.raises(GeometryMissingPhaseAssignmentError):
        CIPHERGeometry(materials=materials, **get_boiler_plate_geometry_args())


def test_geometry_raise_on_missing_phases_second_material_with_first_material_phase_types():
    materials = [
        MaterialDefinition(
            name="mat1",
            properties={},
            phase_types=[
                PhaseTypeDefinition(type_label="1", phases=[0, 1]),
                PhaseTypeDefinition(type_label="2", phases=[2, 3]),
            ],
        ),
        MaterialDefinition(
            name="mat2",
            properties={},
        ),
    ]
    with pytest.raises(GeometryMissingPhaseAssignmentError):
        CIPHERGeometry(materials=materials, **get_boiler_plate_geometry_args())


def test_geometry_raise_on_non_consecutive_voxel_phase():
    with pytest.raises(GeometryVoxelPhaseError):
        CIPHERGeometry(
            voxel_phase=np.array(
                [
                    [0, 0],
                    [2, 2],
                ]
            ),
            size=[1, 1],
            materials=[MaterialDefinition(name="mat1", properties={})],
            interfaces=[],
        )


def test_geometry_single_phase_type_with_specified_phase():
    materials = [MaterialDefinition(name="mat1", properties={}, phases=[0, 1, 2])]
    geom = CIPHERGeometry(
        materials=materials,
        **get_boiler_plate_geometry_args(
            num_phases=3,
            interfaces=[
                InterfaceDefinition(
                    properties={},
                    materials=("mat1", "mat1"),
                    phase_pairs=[[0, 1], [0, 2], [1, 2]],
                )
            ],
        ),
    )
    assert np.all(geom.materials[0].phase_types[0].phases == [0, 1, 2])


def test_geometry_material_phases_concatenation():
    mat = MaterialDefinition(
        name="mat1",
        properties={},
        phase_types=[
            PhaseTypeDefinition(type_label="1", phases=[0, 1, 2]),
            PhaseTypeDefinition(type_label="2", phases=[4, 5]),
        ],
    )
    assert np.all(mat.phases == np.array([0, 1, 2, 4, 5]))


def test_material_raise_on_two_phase_types_but_one_phases():
    with pytest.raises(MaterialPhaseTypePhasesMissingError):
        MaterialDefinition(
            name="mat1",
            properties={},
            phase_types=[
                PhaseTypeDefinition(type_label="1", phases=[1, 3, 2]),
                PhaseTypeDefinition(type_label="2", target_type_fraction=0.6),
            ],
        )


def test_geometry_material_volume_fractions_equal():
    materials = [
        MaterialDefinition(name="mat1", properties={}),
        MaterialDefinition(name="mat2", properties={}),
    ]
    geom = CIPHERGeometry(
        materials=materials,
        **get_boiler_plate_geometry_args(
            interfaces=[
                InterfaceDefinition(materials=("mat1", "mat1"), properties={}),
                InterfaceDefinition(materials=("mat1", "mat2"), properties={}),
                InterfaceDefinition(materials=("mat2", "mat2"), properties={}),
            ]
        ),
    )

    assert [i.target_volume_fraction for i in geom.materials] == [0.5, 0.5]


def test_geometry_material_volume_fractions_expected_remaining():
    materials = [
        MaterialDefinition(name="mat1", properties={}, target_volume_fraction=0.3),
        MaterialDefinition(name="mat2", properties={}),
    ]
    geom = CIPHERGeometry(
        materials=materials,
        **get_boiler_plate_geometry_args(
            interfaces=[
                InterfaceDefinition(materials=("mat1", "mat1"), properties={}),
                InterfaceDefinition(materials=("mat1", "mat2"), properties={}),
                InterfaceDefinition(materials=("mat2", "mat2"), properties={}),
            ]
        ),
    )

    assert geom.materials[1].target_volume_fraction == 0.7


def test_geometry_material_volume_fractions_expected_remaining_two():
    materials = [
        MaterialDefinition(name="mat1", properties={}, target_volume_fraction=0.3),
        MaterialDefinition(name="mat2", properties={}),
        MaterialDefinition(name="mat3", properties={}),
    ]
    geom = CIPHERGeometry(
        materials=materials,
        **get_boiler_plate_geometry_args(
            interfaces=[
                InterfaceDefinition(materials=("mat1", "mat1"), properties={}),
                InterfaceDefinition(materials=("mat1", "mat2"), properties={}),
                InterfaceDefinition(materials=("mat1", "mat3"), properties={}),
                InterfaceDefinition(materials=("mat2", "mat2"), properties={}),
                InterfaceDefinition(materials=("mat2", "mat3"), properties={}),
                InterfaceDefinition(materials=("mat3", "mat3"), properties={}),
            ]
        ),
    )

    assert (
        geom.materials[1].target_volume_fraction
        == geom.materials[2].target_volume_fraction
        == 0.35
    )


def test_geometry_raise_on_duplicate_material_name():
    materials = [
        MaterialDefinition(name="mat1", properties={}),
        MaterialDefinition(name="mat1", properties={}),
    ]
    with pytest.raises(GeometryDuplicateMaterialNameError):
        CIPHERGeometry(materials=materials, **get_boiler_plate_geometry_args())


def test_geometry_raise_on_non_unit_target_volume_fraction():
    materials = [
        MaterialDefinition(
            name="mat1",
            properties={},
            target_volume_fraction=0.9,
        ),
        MaterialDefinition(name="mat2", properties={}, target_volume_fraction=0.3),
    ]
    with pytest.raises(GeometryNonUnitTargetVolumeFractionError):
        CIPHERGeometry(materials=materials, **get_boiler_plate_geometry_args())


def test_geometry_raise_on_excess_target_volume_fraction():
    materials = [
        MaterialDefinition(
            name="mat1",
            properties={},
            target_volume_fraction=1.0,
        ),
        MaterialDefinition(
            name="mat2",
            properties={},
        ),
    ]
    with pytest.raises(GeometryExcessTargetVolumeFractionError):
        CIPHERGeometry(materials=materials, **get_boiler_plate_geometry_args())


def test_geometry_raise_on_excess_target_volume_fraction_all_specified():

    materials = [
        MaterialDefinition(
            name="mat1",
            properties={},
            target_volume_fraction=1.0,
        ),
        MaterialDefinition(
            name="mat2",
            properties={},
            target_volume_fraction=0.2,
        ),
    ]
    with pytest.raises(GeometryNonUnitTargetVolumeFractionError):
        CIPHERGeometry(materials=materials, **get_boiler_plate_geometry_args())


def test_write_input_YAML_round_trip(tmp_path):

    solution_params = {
        "abstol": 0.0001,
        "amrinterval": 25,
        "initblocksize": [1, 1, 1],
        "initcoarsen": 6,
        "initrefine": 7,
        "interpolation": "cubic",
        "maxnrefine": 7,
        "minnrefine": 0,
        "outfile": "out",
        "outputfreq": 100,
        "petscoptions": "-ts_adapt_monitor -ts_rk_type 2a",
        "random_seed": 1579993586,
        "reltol": 0.0001,
        "time": 100000000,
    }

    mat_props = {
        "chemicalenergy": "none",
        "molarvolume": 1e-5,
        "temperature0": 500.0,
    }

    int_props_1 = {
        "energy": {"e0": 5e8},
        "mobility": {"m0": 1e-11},
        "width": 4.0,
    }

    materials = [
        MaterialDefinition(
            name="mat1",
            properties=mat_props,
        ),
        MaterialDefinition(
            name="mat2",
            properties=mat_props,
        ),
    ]

    # Define the interfaces:
    interfaces = [
        InterfaceDefinition(
            materials=("mat1", "mat2"),
            properties=int_props_1,
        ),
        InterfaceDefinition(
            materials=("mat1", "mat1"),
            properties=int_props_1,
        ),
        InterfaceDefinition(
            materials=("mat2", "mat2"),
            properties=int_props_1,
        ),
    ]

    inp = CIPHERInput.from_random_voronoi(
        materials=materials,
        num_phases=500,
        grid_size=[128, 128, 128],
        size=[128, 128, 128],
        components=["ti"],
        outputs=["phaseid", "matid", "interfaceid"],
        solution_parameters=solution_params,
        interfaces=interfaces,
    )

    test_input_path = tmp_path / "round_trip_test.yaml"
    inp.write_yaml(test_input_path)
    inp_reload = CIPHERInput.from_input_YAML_file(test_input_path)
    assert inp == inp_reload
