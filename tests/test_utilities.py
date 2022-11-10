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


def test_get_subset_indices_brute_force():
    for i in range(1, 500):
        for j in range(1, 500):
            out_ij = get_subset_indices(i, j)
            try:
                # check unique indices
                assert len(set(out_ij)) == len(out_ij)
            except AssertionError:
                raise AssertionError(f"{i} grab {j} not unique: {out_ij}.")
            try:
                # check requested length
                if i < j:
                    assert len(out_ij) == i
                else:
                    assert len(out_ij) == j
            except AssertionError:
                raise AssertionError(f"{i} grab {j} length is {len(out_ij)}.")
