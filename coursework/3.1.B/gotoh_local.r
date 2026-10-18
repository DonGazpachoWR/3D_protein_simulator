local_gotoh = function(seq1, seq2, gap, gap_extension, mismatch, match){
        
        build_matrix = function(seq1, seq2){
                matrix_name = matrix(0, nrow = len1 + 1, ncol = len2 + 1) 
                rownames(matrix_name) = c("-", seq1)
                colnames(matrix_name) = c("-", seq2)
                
                return(matrix_name)
        }
        
        
        # Stop conditions
        stopifnot(gap <= 0) # check if penalty negative
        stopifnot(mismatch <= 0)  # check if penalty negative
        stopifnot(match >= 0)  # check if score positive
        stopifnot(gap_extension <= 0) # Added. check if penalty negative
        
        # Initialize col and rownames for matrices
        len1 = nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
        seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
        seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector
        
        # Initialize matrix G for match scores 
        score_G = build_matrix(seq1, seq2)
        # Matrix for horizontal gaps
        score_E = build_matrix(seq1, seq2)
        # Matrix for vertical gaps
        score_F = build_matrix(seq1, seq2)
        # Matrix for maximum scores
        score_V = build_matrix(seq1, seq2)
        
        # Initialize matrix for directions
        dir_G = build_matrix(seq1, seq2)
        dir_E = build_matrix(seq1, seq2)
        dir_F = build_matrix(seq1, seq2)
        
        type = c("dia", "hor", "ver") # Lookup vector
        
        # Build direction matrix
        dir_G[1, 2:ncol(dir_G)] = "hor" 
        dir_G[2:nrow(dir_G), 1] = "ver" 
        dir_E[1, 2:ncol(dir_E)] = "hor"
        dir_F[2:nrow(dir_F), 1] = "ver"
        
        # Compute scores and save moves
        for (i in 2:(len1 + 1)){# for every (initially zero) row
                for (j in 2:(len2 + 1)){# for every (initially zero) col
                        # Score for matrix G (diagonal)
                        score_G[i, j] = ifelse(rownames(score_G)[i] == colnames(score_G)[j], 
                                               score_V[i - 1, j - 1] + match,
                                               score_V[i - 1, j - 1] + mismatch)
                        # Score for matrix E (horizontal)
                        score_E[i, j] = max(score_E[i, j - 1] + gap_extension, 
                                            score_V[i, j - 1] + gap)
                        # Score for matrix F (vertical)
                        score_F[i, j] = max(score_F[i - 1, j] + gap_extension, 
                                            score_V[i - 1, j] + gap)
                        # Score for matrix V (which keep maximum score of the other 3 matrix)
                        score_V[i, j] = max(0, score_G[i, j], score_E[i, j], score_F[i, j])
                        
                        # Direction for matrix G
                        dia = score_G[i, j]
                        hor = score_E[i, j]
                        ver = score_F[i, j]
                        dir_G[i, j] = type[which.max(c(dia, hor, ver))]
                        
                        # Direction for matrix E
                        dir_E[i, j] = ifelse(score_E[i, j - 1] + gap_extension > score_V[i, j - 1] + gap, 
                                             "hor", "ver")
                        
                        # Direction for matrix F
                        dir_F[i, j] = ifelse(score_F[i - 1, j] + gap_extension > score_V[i - 1, j] + gap, 
                                             "ver", "hor")
                        
                }
        } 
        
        
        # Backtracing
        location = as.vector(which(score_V == max(score_V), arr.ind = TRUE))
        i = location[1]
        j = location[2]
        score = score_V[i, j]
        highest_score = score
        align1 = c(); align2 = c() # Note: length of final alignments is unknown at this point
        m_position = "G"
        while(score != 0){
                
                if (m_position == "G"){
                        if(dir_G[i, j] == "dia"){
                                align1 = c(rownames(score_V)[i], align1)
                                align2 = c(colnames(score_V)[j], align2)
                                j = j - 1; i = i - 1  # update indices
                        }
                        else if (dir_G[i, j] == "hor") {
                                m_position = "E" }
                        
                        else if (dir_G[i, j] == "ver") {
                                m_position = "F"} 
                        
                } else if (m_position == "E"){
                        
                        if (dir_E[i, j] == "ver") {
                                m_position = "G"
                        }
                        align1 = c("-", align1) # horizontal movement = gap for seq1
                        align2 = c(colnames(score_V)[j], align2) 
                        j = j - 1 # update indices}
                        
                        
                } else if (m_position == "F"){
                        
                        if (dir_F[i, j] == "hor") {
                                m_position = "G" }
                        align1 = c(rownames(score_V)[i], align1)
                        align2 = c("-", align2) # vertical movement = gap for seq2
                        i = i - 1 # update indices
                }
                
                score = score_V[i, j]
                
                
        }
 
        # Prepare output
        return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2), 
                    score = highest_score, score_V = score_V, score_F, score_E, score_G, dir_G, dir_E, dir_F))}

solution_local_gotoh <- local_gotoh("CTCGTTTCAGAAC", "CGACTCGTTAGAT", gap = -3, gap_extension = -2, mismatch = -4, match = 5)
cat("\nAligned sequences\n")
aligned_seqs_g <- solution_local_gotoh$aligned_seqs
aligned_seqs_g

cat("Score")
score_g <- solution_local_gotoh$score
score_g
score_g_V <- solution_local_gotoh$score_V
write.csv(solution_local_gotoh$score_V, file = "score_matrix_gotoh_local.csv")
