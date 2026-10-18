# This algoritm searches for the best local alignment
smith_waterman = function(seq1, seq2, gap, mismatch, match){
        
        # Stop conditions
        stopifnot(gap <= 0) # check if penalty negative
        stopifnot(mismatch <= 0)  # check if penalty negative
        stopifnot(match >= 0)  # check if score positive
        
        # Initialize col and rownames for matrices
        len1 = nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
        seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
        seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector
        
        # Initialize matrix M (for scores)
        M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
        
        rownames(M) = c("-", seq1) # assign seq chars to matrix names
        colnames(M) = c("-", seq2) # assign seq chars to matrix names
        
        # No gap penalties in first row and column
        
        
        # Initialize matrix D (for directions)
        D = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
        
        rownames(D) = c("-", seq1) # assign seq chars to matrix names
        colnames(D) = c("-", seq2) # assign seq chars to matrix names
        
        D[1, ] = rep("hor") # Fill 1st row with "hor" for horizontal moves
        D[, 1] = rep("ver") # Fill 1st col with "ver" for vertical moves
        
        type = c("dia", "hor", "ver") # Lookup vector
        
        # Compute scores and save moves
        for (i in 2:(len1 + 1)){# for every (initially zero) row
                for (j in 2:(len2 + 1)){# for every (initially zero) col
                        hor = M[i, j - 1] + gap # horizontal move = gap for seq1
                        ver = M[i - 1, j] + gap # vertical move = gap for seq2
                        dia = ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch) 
                                     M[i - 1, j - 1] + match, 
                                     M[i - 1, j - 1] + mismatch)
                        # We add 0 for choosing the max, so we keep the value 0 as the lowest score in the matrix
                        M[i, j] = max(0, dia, hor, ver) # Save current (best) score in M
                        D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
                }
        } 
        
        
        # Backtracing
        align1 = c(); align2 = c() # Note: length of final alignments is unknown at this point
        
        # The backtracking starts in the maximum score, not the last. It ends when it finds a 0
        location = as.vector(which(M == max(M), arr.ind = TRUE))
        i = location[1]
        j = location[2]
        score = M[i, j]
        highest_score = score
        
        while(score != 0){
                
                if(D[i, j] == "dia") {
                        align1 = c(rownames(M)[i], align1)
                        align2 = c(colnames(M)[j], align2)
                        j = j - 1; i = i - 1  # update indices
                } else if (D[i, j] == "ver") {
                        align1 = c(rownames(M)[i], align1)
                        align2 = c("-", align2) # vertical movement = gap for seq2
                        i = i - 1 # update indices
                } else if (D[i, j] == "hor") {
                        align1 = c("-", align1) # horizontal movement = gap for seq1
                        align2 = c(colnames(M)[j], align2) 
                        j = j - 1 # update indices
                } 
                score = M[i, j]
        }
        
        # Prepare output
        return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2),
                    score = highest_score, score_matrix = M, movement_matrix = D))
        
}

solution_s <- smith_waterman("CTCGTTTCAGAAC", "CGACTCGTTAGAT", gap = -3, mismatch = -4, match = 5)
cat("\nAligned sequences\n")
aligned_seqs_s <- solution_s$aligned_seqs
aligned_seqs_s
cat("\nMovement matrix\n")
movement_matrix_s <- solution_s$movement_matrix
movement_matrix_s
cat("\nScore matrix\n")
score_matrix_s <- solution_s$score_matrix
score_matrix_s
cat("Score")
score_s <- solution_s$score
score_s
write.csv(score_matrix_s, file = "score_matrix_Smith_waterman.csv")