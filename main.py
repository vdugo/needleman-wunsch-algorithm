def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
	"""
    Perform global alignment of two DNA sequences using the Needleman-Wunsch algorithm.

    Args:
        seq1: A string representing the first DNA sequence to align.
        seq2: A string representing the second DNA sequence to align.
        match_score: The score to assign to matches between two nucleotides (default 1).
        mismatch_score: The score to assign to mismatches between two nucleotides (default -1).
        gap_penalty: The score to assign for inserting a gap in the alignment (default -1).

    Returns:
        A tuple containing the two aligned DNA sequences and their alignment score.

    The function constructs a scoring matrix by computing the optimal alignment score for each pair of
    nucleotides in the two sequences. It then traces back from the bottom-right corner of the matrix to
    construct the optimal alignment of the two sequences.

    This implementation assumes that the input sequences consist only of the characters 'A', 'C', 'G', 'T',
    and will raise an error if this is not the case.
    """
	
	# Initialize scoring matrix with size (n + 1) x (m + 1)
	n, m = len(seq1), len(seq2)
	M = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

	# Initialize the first row and column with gap penalties
	for i in range(n + 1):
		M[i][0] = i * gap_penalty
	for j in range(m + 1):
		M[0][j] = j * gap_penalty

	# Filling the matrix
	for i in range(1, n + 1):
		for j in range(1, m + 1):
			# if the nucleotides are aligned between the two sequences
			if seq1[i - 1] == seq2[j - 1]:
			    diag_score = match_score
			else:
			    diag_score = mismatch_score
			
			if diag_score < 0:
			    score_diag = M[i - 1][j - 1] - abs(diag_score)
			else:
			    score_diag = M[i - 1][j - 1] + diag_score

			# calculate the scores for each possible move
			score_up = M[i - 1][j] + gap_penalty
			score_down = M[i][j - 1] + gap_penalty

			# choose the maximum score to fill the current cell
			M[i][j] = max(score_diag, score_up, score_down)

	# traceback
	aligned_seq1 = []
	aligned_seq2 = []
	# start from the bottom right of the matrix
	i, j = n, m

	while i > 0 or j > 0:
		# if the move is diagonal
		if i > 0 and j > 0 and M[i][j] == M[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
			aligned_seq1.append(seq1[i-1])
			aligned_seq2.append(seq2[j-1])
			i -= 1
			j -= 1
		# else if the move is vertical
		elif i > 0 and M[i][j] == M[i-1][j] + gap_penalty:
			aligned_seq1.append(seq1[i-1])
			aligned_seq2.append('-')
			i -= 1
		# else the move is horizontal
		elif j > 0 and M[i][j] == M[i][j-1] + gap_penalty:
			aligned_seq1.append('-')
			aligned_seq2.append(seq2[j-1])
			j -= 1


	# reverse the aligned sequences
	aligned_seq1 = ''.join(aligned_seq1[::-1])
	aligned_seq2 = ''.join(aligned_seq2[::-1])

	return aligned_seq1, aligned_seq2, M[n][m]

seq1, seq2, score = needleman_wunsch("GATTACA", "GCATGCU")

print(seq1)
print(seq2)
print(score)
			

	
