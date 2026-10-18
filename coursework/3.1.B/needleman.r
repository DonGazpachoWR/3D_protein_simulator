# New version of https://gist.github.com/juliuskittler/ed53696ac1e590b413aac2dddf0457f6 
# with less penalty for extending the gap

needleman = function(seq1, seq2, gap, mismatch, match){
        
        # Stop conditions
        stopifnot(gap <= 0) # check if penalty negative
        stopifnot(mismatch <= 0)  # check if penalty negative
        # Initialize col and rownames for matrices
        len1 = nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
        seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
        seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector
        
        # Initialize matrix M (for scores)
        M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
        
        rownames(M) = c("-", seq1) # assign seq chars to matrix names
        colnames(M) = c("-", seq2) # assign seq chars to matrix names
        
        M[1, ] = cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
        M[, 1] = cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites
        
        
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
                        M[i, j] = max(dia, hor, ver) # Save current (best) score in M
                        D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
                }
        } 
       
        
        # Backtracing
        align1 = c(); align2 = c() # Note: length of final alignments is unknown at this point
        
        while(i > 1 || j > 1){
                
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
                
        }

        # Prepare output
        return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2),
                    score = M[nrow(M), ncol(M)], score_matrix = M, movement_matrix = D))
        
}

# Needleman algorithm
solution_n <- needleman("CTCGTTTCAGAAC", "CGACTCGTTAGAT", gap = -3, mismatch = -4, match = 5)
cat("\nAligned sequences\n")
aligned_seqs_n <- solution_n$aligned_seqs
aligned_seqs_n
cat("\nMovement matrix\n")
movement_matrix_n <- solution_n$movement_matrix
movement_matrix_n
cat("\nScore matrix\n")
score_matrix_n <- solution_n$score_matrix
score_matrix_n
cat("Score")
score_n <- solution_n$score
score_n
write.csv(score_matrix_n, file = "score_matrix_needleman.csv")






