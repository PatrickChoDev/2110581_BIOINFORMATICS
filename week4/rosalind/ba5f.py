import numpy as np


def local_alignment_score(s1, s2):
    alphabet_index = [c for c in "ACDEFGHIKLMNPQRSTVWY"]
    # Scoring matrix (PAM250)
    scoring_matrix = np.array(
        [
            [2, -2, 0, 0, -3, 1, -1, -1, -1, -2, -1, 0, 1, 0, -2, 1, 1, 0, -6, -3],
            [
                -2,
                12,
                -5,
                -5,
                -4,
                -3,
                -3,
                -2,
                -5,
                -6,
                -5,
                -4,
                -3,
                -5,
                -4,
                0,
                -2,
                -2,
                -8,
                0,
            ],
            [0, -5, 4, 3, -6, 1, 1, -2, 0, -4, -3, 2, -1, 2, -1, 0, 0, -2, -7, -4],
            [0, -5, 3, 4, -5, 0, 1, -2, 0, -3, -2, 1, -1, 2, -1, 0, 0, -2, -7, -4],
            [-3, -4, -6, -5, 9, -5, -2, 1, -5, 2, 0, -3, -5, -5, -4, -3, -3, -1, 0, 7],
            [1, -3, 1, 0, -5, 5, -2, -3, -2, -4, -3, 0, 0, -1, -3, 1, 0, -1, -7, -5],
            [-1, -3, 1, 1, -2, -2, 6, -2, 0, -2, -2, 2, 0, 3, 2, -1, -1, -2, -3, 0],
            [-1, -2, -2, -2, 1, -3, -2, 5, -2, 2, 2, -2, -2, -2, -2, -1, 0, 4, -5, -1],
            [-1, -5, 0, 0, -5, -2, 0, -2, 5, -3, 0, 1, -1, 1, 3, 0, 0, -2, -3, -4],
            [-2, -6, -4, -3, 2, -4, -2, 2, -3, 6, 4, -3, -3, -2, -3, -3, -2, 2, -2, -1],
            [-1, -5, -3, -2, 0, -3, -2, 2, 0, 4, 6, -2, -2, -1, 0, -2, -1, 2, -4, -2],
            [0, -4, 2, 1, -3, 0, 2, -2, 1, -3, -2, 2, 0, 1, 0, 1, 0, -2, -4, -2],
            [1, -3, -1, -1, -5, 0, 0, -2, -1, -3, -2, 0, 6, 0, 0, 1, 0, -1, -6, -5],
            [0, -5, 2, 2, -5, -1, 3, -2, 1, -2, -1, 1, 0, 4, 1, -1, -1, -2, -5, -4],
            [-2, -4, -1, -1, -4, -3, 2, -2, 3, -3, 0, 0, 0, 1, 6, 0, -1, -2, 2, -4],
            [1, 0, 0, 0, -3, 1, -1, -1, 0, -3, -2, 1, 1, -1, 0, 2, 1, -1, -2, -3],
            [1, -2, 0, 0, -3, 0, -1, 0, 0, -2, -1, 0, 0, -1, -1, 1, 3, 0, -5, -3],
            [0, -2, -2, -2, -1, -1, -2, 4, -2, 2, 2, -2, -1, -2, -2, -1, 0, 4, -6, -2],
            [
                -6,
                -8,
                -7,
                -7,
                0,
                -7,
                -3,
                -5,
                -3,
                -2,
                -4,
                -4,
                -6,
                -5,
                2,
                -2,
                -5,
                -6,
                17,
                0,
            ],
            [
                -3,
                0,
                -4,
                -4,
                7,
                -5,
                0,
                -1,
                -4,
                -1,
                -2,
                -2,
                -5,
                -4,
                -4,
                -3,
                -3,
                -2,
                0,
                10,
            ],
        ]
    )

    # Indel penalty
    indel_penalty = 5

    # Initialize the score matrix
    score_matrix = np.zeros((len(s1) + 1, len(s2) + 1))

    # Initialize variables to store the maximum score and its indices
    max_score = 0
    max_i = 0
    max_j = 0

    # Fill the score matrix
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):

            match_score = (
                score_matrix[i - 1, j - 1]
                + scoring_matrix[
                    alphabet_index.index(s1[i - 1]),
                    alphabet_index.index(s2[j - 1]),
                ]
            )
            delete_score = score_matrix[i - 1, j] - indel_penalty
            insert_score = score_matrix[i, j - 1] - indel_penalty
            score_matrix[i, j] = max(0, match_score, delete_score, insert_score)

            # Update the maximum score and its indices
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_i = i
                max_j = j
    print(score_matrix)
    # Traceback to construct the local alignment
    alignment_s1 = ""
    alignment_s2 = ""
    i, j = max_i, max_j
    while score_matrix[i, j] > 0:
        if (
            score_matrix[i, j]
            == score_matrix[i - 1, j - 1]
            + scoring_matrix[
                alphabet_index.index(s1[i - 1]), alphabet_index.index(s2[j - 1])
            ]
        ):
            alignment_s1 = s1[i - 1] + alignment_s1
            alignment_s2 = s2[j - 1] + alignment_s2
            i -= 1
            j -= 1
        elif score_matrix[i, j] == score_matrix[i - 1, j] - indel_penalty:
            alignment_s1 = s1[i - 1] + alignment_s1
            alignment_s2 = "-" + alignment_s2
            i -= 1
        else:
            alignment_s1 = "-" + alignment_s1
            alignment_s2 = s2[j - 1] + alignment_s2
            j -= 1

    return max_score, alignment_s1, alignment_s2


with open("data/rosalind_ba5f.txt", "r") as f:
    s1 = f.readline().strip()
    s2 = f.readline().strip()

max_score, alignment_s1, alignment_s2 = local_alignment_score(s1, s2)
with open("rosalind_ba5f.out", "w") as f:
    f.write(str(int(max_score)) + "\n")
    f.write(alignment_s1 + "\n")
    f.write(alignment_s2 + "\n")

print(int(max_score))
print(alignment_s1)
print(alignment_s2)
