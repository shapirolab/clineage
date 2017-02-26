ms_28727_a_id = 1
ms_28727_b_id = 2
ms_28734_a_id = 3
ms_adj_ms_1_a = 4
ms_adj_ms_2_a = 5
ms_adj_ms_2_b = 6
amplicon_28727_id = 1
amplicon_28734_id = 2
pu_adj_ms_1 = 3
pu_adj_ms_2 = 4


def repeat_number(i):  # Just to give the int some context
    return i


MS_HISTOGRAMS_DICT = {
    amplicon_28734_id: {
        ms_28734_a_id: {
            frozenset([(1., 17)]): {
                repeat_number(7): 1,
                repeat_number(8): 6,
                repeat_number(9): 8,
                repeat_number(10): 18,
                repeat_number(11): 24,
                repeat_number(12): 117,
                repeat_number(13): 329,
                repeat_number(14): 824,
                repeat_number(15): 2081,
                repeat_number(16): 3899,
                repeat_number(17): 4071,
                repeat_number(18): 346,
                repeat_number(19): 30,
                repeat_number(20): 5,
                repeat_number(21): 1,
            },
            frozenset([(0.6, 17), (0.4, 27)]): {  # biallelic example
                repeat_number(6): 1,
                repeat_number(7): 1,
                repeat_number(8): 3,
                repeat_number(9): 5,
                repeat_number(10): 6,
                repeat_number(11): 12,
                repeat_number(12): 27,
                repeat_number(13): 83,
                repeat_number(14): 189,
                repeat_number(15): 398,
                repeat_number(16): 761,
                repeat_number(17): 792,
                repeat_number(18): 77,
                repeat_number(19): 24,
                repeat_number(20): 28,
                repeat_number(21): 62,
                repeat_number(22): 113,
                repeat_number(23): 173,
                repeat_number(24): 311,
                repeat_number(25): 394,
                repeat_number(26): 375,
                repeat_number(27): 220,
                repeat_number(28): 47,
                repeat_number(29): 5,
            },
        },
    },
}
