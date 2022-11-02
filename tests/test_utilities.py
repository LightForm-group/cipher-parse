from cipher_parse.utilities import get_subset_indices


def test_get_subset_indices():
    assert get_subset_indices(7, 0) == []
    assert get_subset_indices(7, 1) == [0]
    assert get_subset_indices(7, 2) == [0, 6]
    assert get_subset_indices(7, 3) == [0, 3, 6]
    assert get_subset_indices(7, 4) == [0, 2, 4, 6]
    assert get_subset_indices(7, 5) == [0, 2, 3, 5, 6]
    assert get_subset_indices(7, 6) == [0, 1, 2, 3, 4, 6]
    assert get_subset_indices(7, 7) == [0, 1, 2, 3, 4, 5, 6]
    assert get_subset_indices(7, 8) == [0, 1, 2, 3, 4, 5, 6]

    assert get_subset_indices(9, 2) == [0, 8]
    assert get_subset_indices(9, 3) == [0, 4, 8]
    assert get_subset_indices(9, 4) == [0, 3, 6, 8]
    assert get_subset_indices(9, 5) == [0, 2, 4, 6, 8]
    assert get_subset_indices(9, 6) == [0, 2, 3, 5, 6, 8]
    assert get_subset_indices(9, 7) == [0, 1, 2, 4, 5, 6, 8]
    assert get_subset_indices(9, 8) == [0, 1, 2, 3, 4, 5, 6, 8]
