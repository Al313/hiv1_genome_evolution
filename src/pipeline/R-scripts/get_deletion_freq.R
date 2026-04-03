


deletion_junctions <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/misc/nef_deletions_list.tsv",
                                sep = "",
                                header = TRUE)


row <- 13
psg <- 620

df <- data.frame()
line_conversion <- c("MT-2_1" = "13", "MT-2_2" = "14", "MT-4_1" = "15", "MT-4_2" = "16")
for (row in 7:7){
    print(row)
    start_psg <- deletion_junctions$start_psg[row]
    end_psg <- deletion_junctions$end_psg[row]
    exp_line <- deletion_junctions$exp_line[row]
    exp_line_no <- as.character(line_conversion[exp_line])


    start_pos <- deletion_junctions$start_pos[row] + deletion_junctions$offset_length[row]
    start_pos_left <- deletion_junctions$start_pos[row] - 1 - deletion_junctions$offset_length[row]
    end_pos <- deletion_junctions$end_pos[row] - deletion_junctions$offset_length[row]
    end_pos_right <- deletion_junctions$end_pos[row] + 1 + deletion_junctions$offset_length[row]

    for (psg in seq(start_psg, end_psg, 10)){
        print(psg)
        file <- list.files(
            path = paste0("/Users/alimos313/Desktop/scrap/bam-portal/", exp_line_no),
            pattern = paste0("VP", psg, ".*\\.bam$"),
            full.names = TRUE
        )

        cov_all <- c()
        for (pos in c(start_pos, start_pos_left, end_pos, end_pos_right)){
            depth <- system2(
                "samtools",
                args = c(
                    "depth",
                    "-q", "35",
                    "-d", "1000000",
                    "-r", paste0("NL43_ann_wk0virusPassRef_plasmid:", pos, "-", pos),
                    file
                ),
                stdout = TRUE
            )

            if (length(depth) == 0){
                cov <- 0
            } else {
                cov <- str_split(depth, pattern = "\t")[[1]][3]
            }
            cov_all <- c(cov_all, cov)
        }
        if (is.na(deletion_junctions$adjacent_deletion[row])){
            NULL
        } else if (deletion_junctions$adjacent_deletion[row] == "start"){
            cov_all[1:2] <- NA
        } else if (deletion_junctions$adjacent_deletion[row] == "end"){
            cov_all[3:4] <- NA
        }
        
        cov_all <- c(as.character(deletion_junctions[row,c(1:2)]), psg ,as.character(deletion_junctions[row,c(3:8)]), cov_all)
        df <- rbind(df, cov_all)
    }
}

colnames(df) <- c("exp_line", "deletion_id", "passage", "start_psg", "end_psg", "start_pos", "end_pos", "offset_length", 
    "adjacent_deletion", "start_pos_cov", "start_pos_left_cov", "end_pos_cov", "end_pos_right_cov")

saveRDS(df, file = paste0(wd,"results/tables/misc/nef_deletions_list_processed.Rds"))


