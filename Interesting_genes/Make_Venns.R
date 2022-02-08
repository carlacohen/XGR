#Aim to make Venn diagrams of interesting genes
library(VennDiagram)

### Read in data

file_list <- list.files(".", pattern="*.txt")

files <- list()

#Use scan to read in text files as a vector

for(i in 1:length(file_list)){
    files[[i]] <- scan(file_list[[i]], what = "character")
}
head(files [[1]])

#Name the files based on the list
names(files) <- file_list

#Read in CD14-related files only

# Make list of CD14 files

file_list.CD14 <- list.files(".", pattern = "CD14")
file_list.CD14

# Read in CD14 files as vector

files.CD14 <- list()

for(i in 1:length(file_list.CD14)){
    files.CD14[[i]] <- scan(file_list.CD14[[i]], what = "character")
    }
head(files.CD14 [[1]])
file_list.CD14

#Remove NAs from lists
test_list <- lapply(files.CD14, na.omit)

#Generate Venn diagram

venn.diagram(
    x = files.CD14, 
    category.names = test_list,
    filename = "Venn.CD14.png",
    output= TRUE
)
