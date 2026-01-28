# Homework 2: Basic Bash and R
## Questions
1. Navigate to your home directory and create a directory called homework-2 then navigate into it. Then create 2 directories labeled scripts and data. In scripts create a file called hello. Once you do that navigate out back to the homework-2 directory and delete the data directory.
2. In Rstudio create a file called scoreboard.R. In that file paste the following:
   ```ruby
   stats_matrix <- matrix(
        c(11, 2, 9,
        30, 5, 0),
        nrow = 2,
        byrow = TRUE
      )
    
   colnames(stats_matrix) <- c("points", "rebounds", "assists")
    
   stats_df <- data.frame(
        points = c(11, 30),
        rebounds = c(2, 5),
        assists = c(9, 0)
      )
   
   ```
    Now once that is pasted I want you to answer the following questions:
    - Print the first players points
    - Print the points for both players.
    - What is the difference between:
        - stats_matrix[, "rebounds"]
        - stats_df[, "rebounds"]
        - stats_df["rebounds"]
        - stats_df$rebounds
        - stats_df[["rebounds"]]

3. Take the final file and make a copy of it in the homework-2/scripts folder and add #!/usr/bin/env Rscript to the top of the file to make it executable. We want to make the file readable and executable to run as well as all parent directories.
   - What 3-digit chmod code should you use with chmod to achieve the permissions described above?
   - What does each component of the code mean?
  


## Answer Key
1. Command list:
    ```ruby
    (base) gautham@dhcp-172-31-162-096 Fundamentals_of_Informatics % mkdir homework-2
    (base) gautham@dhcp-172-31-162-096 Fundamentals_of_Informatics % cd homework-2 
    (base) gautham@dhcp-172-31-162-096 homework-2 % mkdir scripts
    (base) gautham@dhcp-172-31-162-096 homework-2 % mkdir data
    (base) gautham@dhcp-172-31-162-096 homework-2 % cd scripts 
    (base) gautham@dhcp-172-31-162-096 scripts % touch hello
    (base) gautham@dhcp-172-31-162-096 scripts % ls
    hello
    (base) gautham@dhcp-172-31-162-096 scripts % cd ../  
    (base) gautham@dhcp-172-31-162-096 homework-2 % rmdir data
    (base) gautham@dhcp-172-31-162-096 homework-2 % ls
    scripts
    ```
2. Code pasted in RStudion
   - ```
     > stats_df[1, "points"]
      [1] 11
     ```
   - ```
     > stats_df[, "points"]
      [1] 11 30
     ```
   - ```
     > stats_matrix[, "rebounds"]
      [1] 2 5
     > stats_df[, "rebounds"]
      [1] 2 5
     > stats_df["rebounds"]
       rebounds
      1        2
      2        5
     > stats_df$rebounds
      [1] 2 5
     > stats_df[["rebounds"]]
      [1] 2 5
     ```
     stats_matrix[, "rebounds"] returns a numeric vector as we are extracting from a matrix object. We are selecting one column of the matrix which would be a vector. stats_df[, "rebounds"] also returns a numeric vector because indexing a dataframe would take just the indexed column(s). stats_df["rebounds"] returns a dataframe with only the rebounds column as we specified. stats_df$rebounds also returns a vector as its grabbing the rebounds column from stats_df. stats_df[["rebounds"]] returns a vector and is considered the most explicit way to extract the rebounds column.
3. Code:
   ```ruby
      1 #!/usr/bin/env Rscript
      2 stats_matrix <- matrix(
      3   c(11, 2, 9,
      4     30, 5, 0),
      5   nrow = 2,
      6   byrow = TRUE
      7 )
      8 
      9 colnames(stats_matrix) <- c("points", "rebounds", "assists")
     10 
     11 stats_df <- data.frame(
     12   points = c(11, 30),
     13   rebounds = c(2, 5),
     14   assists = c(9, 0)
     15 )
     16 
     17 print(stats_df[1, "points"])
     18 print(stats_df[, "points"])
   ```
   chmod 755 scoreboard.R
   
   The code 755 means 7 applies to the owner giving them read, write and execute permissions. 5 applies to the group giving them read and execute permissions, but they cannot write. The last 5 is for others also giving them read and execute permissions, but they cannot write.
