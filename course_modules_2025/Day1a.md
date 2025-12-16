# Polygenic Risk Score Analyses Workshop 2025


# Day 1a: GWAS & relevant Statistics


## Introduction to Bash

 Most software in Bioinformatics and Statistical Genetics need to be
 run in a Unix environment (e.g. Linux or Mac OS) and most
 high-performance computer clusters run Unix systems. Therefore,
 although there are alternatives available on Windows (command line,
 Linux subsystems or Virtual Machines), it will be highly beneficial to
 become familiar with performing research in a Unix-only environment.

## What is a shell and why should I care?

A *shell* is a computer program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard/touchscreen combination.

There are many reasons to learn about the shell:

- Many bioinformatics tools can only be used through a command line interface. Many more
  have features and parameter options which are not available in the GUI.
  BLAST is an example. Many of the advanced functions are only accessible
  to users who know how to use a shell.
- The shell makes your work less boring. In bioinformatics you often need to repeat tasks with a large number of files. With the shell, you can automate those repetitive tasks and leave you free to do more exciting things.
- The shell makes your work less error-prone. When humans do the same thing a hundred different times
  (or even ten times), they're likely to make a mistake. Your computer can do the same thing a thousand times
  with no mistakes.
- The shell makes your work more reproducible. When you carry out your work in the command-line
  (rather than a GUI), your computer keeps a record of every step that you've carried out, which you can use
  to re-do your work when you need to. It also gives you a way to communicate unambiguously what you've done,
  so that others can inspect or apply your process to new data.
- Many bioinformatic tasks require large amounts of computing power and can't realistically be run on your
  own machine. These tasks are best performed using remote computers or cloud computing, which can only be accessed
  through a shell.

In this lesson you will learn how to use the command line interface to move around in your file system.

## How to access the shell

On a Mac or Linux machine, you can access a shell through a program called "Terminal", which is already available
on your computer. The Terminal is a window into which we will type commands. If you're using Windows,
you'll need to download a separate program to access the shell.

To save time, we are going to be working on a remote server where all the necessary data and software available.
When we say a 'remote server', we are talking about a computer that is not the one you are working on right now.
You will access the Carpentries remote server where everything is prepared for the lesson.
We will learn the basics of the shell by manipulating some data files. Some of these files are very large
, and would take time to download to your computer.
We will also be using several bioinformatic packages in later lessons and installing all of the software
would take up time even more time. A 'ready-to-go' server lets us focus on learning.

## Moving around the File System

The part of the operating system that manages files and directories
is called the **file system**.
It organizes our data into files,
which hold information,
and directories (also called "folders"),
which hold files or other directories.

Several commands are frequently used to create, inspect, rename, and delete files and directories.
```bash
$
```

The dollar sign is a **prompt**, which shows us that the shell is waiting for input;
your shell may use a different character as a prompt and may add information before
the prompt. When typing commands, either from these lessons or from other sources,
do not type the prompt, only the commands that follow it.

Let's find out where we are by running a command called `pwd`
(which stands for "print working directory").
At any moment, our **current working directory**
is our current default directory,
i.e.,
the directory that the computer assumes we want to run commands in,
unless we explicitly specify something else.
Here,
the computer's response is `/data`,
which is the top level directory within our cloud system:

```bash
$ pwd
```

To begin our practical, please open up a \"terminal\" on your computer
 (on a Mac this is stored in Applications/Utilities/).

 We can change our directory using the following command:

       cd \<Path>\

 where *\<Path>\* is the path to the target directory.


 Some common usage of cd includes

 

       cd ~/ # will bring you to your home directory

       cd ../ # will bring you to the parent directory (up one level)

       cd XXX # will bring you to the XXX directory, so long as it is in the current directory

 As an example, we can move to the **data** directory by
 typing:

       cd data/



## Looking at the Current Directory

 Next we can move into the **~/data/Day1a_Data/Day1a_Data** folder (from the data/ folder type: cd Day1a_Data/Day1a_Data). We can list out

 Let's look at how our file system is organized. We can see what files and subdirectories are in this directory by running `ls`,
which stands for "listing":

```bash
$ ls
```

`ls` prints the names of the files and directories in the current directory in
alphabetical order,
arranged neatly into columns.

We can make the `ls` output more comprehensible by using the **flag** `-F`,
which tells `ls` to add a trailing `/` to the names of directories:

```bash
$ ls -F
```
Anything with a "/" after it is a directory. Things with a "\*" after them are programs. If
there are no decorations, it's a file.

`ls` has lots of other options. To find out what they are, we can type:

```bash
$ man ls
```

`man` (short for manual) displays detailed documentation (also referred as man page or man file)
for `bash` commands. It is a powerful resource to explore `bash` commands, understand
their usage and flags. Some manual files are very long. You can scroll through the
file using your keyboard's down arrow or use the <kbd>Space</kbd> key to go forward one page
and the <kbd>b</kbd> key to go backwards one page. When you are done reading, hit <kbd>q</kbd>
to quit.


to read and write to the file.
 the folder content by typing: 

        ls

 For ls, there are a number of additional Unix command options that you
 can append to it to get additional information, for example:

        ls -l  # shows files as list

        ls -lh  # shows files as a list with human readable format

        ls -lt  # shows the files as a list sorted by time-last-edited

        ls -lS  # shows the files as a list sorted by size

The additional information given includes the name of the owner of the file,
when the file was last modified, and whether the current user has permission

## Full vs. Relative Paths

The `cd` command takes an argument which is a directory
name. Directories can be specified using either a *relative* path or a
full *absolute* path. The directories on the computer are arranged into a
hierarchy. The full path tells you where a directory is in that
hierarchy. Navigate to the home directory, then enter the `pwd`
command.

```bash
$ cd  
$ pwd  
```

You will see:

```output
/home/data
```

This is the full name of your home directory. This tells you that you
are in a directory called `data`, which sits inside a directory called
`home` which sits inside the very top directory in the hierarchy. The
very top of the hierarchy is a directory called `/` which is usually
referred to as the *root directory*. So, to summarize: `data` is a
directory in `home` which is a directory in `/`. More on `root` and
`home` in the next section.

Now enter the following command:

```bash
$ cd /home/Data/Day_1a/
```

This jumps forward multiple levels to the `Day_1a` directory.
Now go back to the home directory.

```bash
$ cd
```

You can also navigate to the `Day_1a` directory using:

```bash
$ cd /Data/Day_1a/
```

These two commands have the same effect, they both take us to the `Day_1a` directory.
The first uses the absolute path, giving the full address from the home directory. The
second uses a relative path, giving only the address from the working directory. A full
path always starts with a `/`. A relative path does not.

A relative path is like getting directions from someone on the street. They tell you to
"go right at the stop sign, and then turn left on Main Street". That works great if
you're standing there together, but not so well if you're trying to tell someone how to
get there from another country. A full path is like GPS coordinates. It tells you exactly
where something is no matter where you are right now.

You can usually use either a full path or a relative path depending on what is most convenient.
If we are in the home directory, it is more convenient to enter the full path.
If we are in the working directory, it is more convenient to enter the relative path
since it involves less typing.

Over time, it will become easier for you to keep a mental note of the
structure of the directories that you are using and how to quickly
navigate amongst them.


## Counting Number of Lines in File

 We can also count the number of lines in a file with the following
 command (where

 *\<file>\* is the file of interest):

        wc -l <file>


 Often we would like to store the output of a command, which we can do
 by *redirecting* the output of the command to a file. For example, we
 can redirect the count of the **GIANT_Height.txt** to **giant_count**
 using the following command:

     wc -l GIANT_Height.txt > giant_count.txt

## Search File Content

 Another common task is to search for specific words or characters in a
 file (e.g. does this file contain our gene of interest?). This can be
 performed using the "grep" command as follows:

       grep <string> file

 For example, to check if the Single Nucleotide Polymorphism
 (SNP) *rs10786427* is present in **GIANT_Height.txt**,
 we can do:

       grep rs10786427 GIANT_Height.txt

 In addition, grep allows us to check if patterns contained in one file
 can be found in another file. For example, if we want to extract a
 subset of samples from the phenotype file (e.g. extract the list of
 samples in **Data/Day_1a/TAR.height**), we can do:

        grep -f Select.sample TAR.height

 An extremely useful feature of the terminal is chaining multiple
 commands into one command, which we call ***piping*** .

 For example, we can use piping to count the number of samples in
 **Select.sample**

 that were found in **TAR.height** in a single command, as follows:

  ```bash 
   grep -f Select.sample TAR.height | wc -l
```

## Creating, moving, copying, and removing

Now we can move around in the file structure, look at files, and search files. But what if we want to copy files or move
them around or get rid of them? Most of the time, you can do these sorts of file manipulations without the command line,
but there will be some cases (like when you're working with a remote computer like we are for this lesson) where it will be
impossible. You'll also find that you may be working with hundreds of files and want to do similar manipulations to all
of those files. In cases like this, it's much faster to do these operations at the command line.

### Copying Files

When working with computational data, it's important to keep a safe copy of that data that can't be accidentally overwritten or deleted.
For this lesson, our raw data is our FASTQ files.  We don't want to accidentally change the original files, so we'll make a copy of them
and change the file permissions so that we can read from, but not write to, the files.

First, let's make a copy of one of our FASTQ files using the `cp` command.

Navigate to the `shell_data/untrimmed_fastq` directory and enter:

```bash
$ cp GIANT_Height.txt GIANT_Height-copy.txt
$ ls -F
```

```output
GIANT_Height.txt GIANT_Height-copy.txt
```

We now have two copies of the `GIANT_Height.txt` file, one of them named `GIANT_Height-copy.txt`. We'll move this file to a new directory
called `backup` where we'll store our backup data files.

### Creating Directories

The `mkdir` command is used to make a directory. Enter `mkdir`
followed by a space, then the directory name you want to create:

```bash
$ mkdir backup
```

### Moving / Renaming

We can now move our backup file to this directory. We can
move files around using the command `mv`:

```bash
$ mv GIANT_Height-copy.txt backup
$ ls backup
```

```output
GIANT_Height-copy.txt
```

The `mv` command is also how you rename files. Let's rename this file to make it clear that this is a backup:

```bash
$ cd backup
$ mv GIANT_Height-copy.txt GIANT_Height-backup.txt
$ ls
```

```output
GIANT_Height-backup.txt
```

### File Permissions

We've now made a backup copy of our file, but just because we have two copies, it doesn't make us safe. We can still accidentally delete or
overwrite both copies. To make sure we can't accidentally mess up this backup file, we're going to change the permissions on the file so
that we're only allowed to read (i.e. view) the file, not write to it (i.e. make new changes).

View the current permissions on a file using the `-l` (long) flag for the `ls` command:

```bash
$ ls -l
```


The first part of the output for the `-l` flag gives you information about the file's current permissions. There are ten slots in the
permissions list. The first character in this list is related to file type, not permissions, so we'll ignore it for now. The next three
characters relate to the permissions that the file owner has, the next three relate to the permissions for group members, and the final
three characters specify what other users outside of your group can do with the file. We're going to concentrate on the three positions
that deal with your permissions (as the file owner).

![](fig/rwx_figure.svg){alt='Permissions breakdown'}

Here the three positions that relate to the file owner are `rw-`. The `r` means that you have permission to read the file, the `w`
indicates that you have permission to write to (i.e. make changes to) the file, and the third position is a `-`, indicating that you
don't have permission to carry out the ability encoded by that space (this is the space where `x` or executable ability is stored, we'll
talk more about this in [a later lesson](05-writing-scripts.md)).

Our goal for now is to change permissions on this file so that you no longer have `w` or write permissions. We can do this using the `chmod` (change mode) command and subtracting (`-`) the write permission `-w`.

```bash
$ chmod -w GIANT_Height-backup.txt
$ ls -l 
```


### Removing

To prove to ourselves that you no longer have the ability to modify this file, try deleting it with the `rm` command:

```bash
$ rm GIANT_Height-backup.txt
```

You'll be asked if you want to override your file permissions:

```output
rm: remove write-protected regular file ‘GIANT_Height-backup.txt'? 
```

You should enter `n` for no. If you enter `n` (for no), the file will not be deleted. If you enter `y`, you will delete the file. This gives us an extra
measure of security, as there is one more step between us and deleting our data files.

**Important**: The `rm` command permanently removes the file. Be careful with this command. It doesn't
just nicely put the files in the Trash. They're really gone.

By default, `rm` will not delete directories. You can tell `rm` to
delete a directory using the `-r` (recursive) option. Let's delete the backup directory
we just made.

Enter the following command:

```bash
$ cd ..
$ rm -r backup
```

This will delete not only the directory, but all files within the directory. If you have write-protected files in the directory,
you will be asked whether you want to override your permission settings.



## Filtering and Reshuﬄing Files

 A very powerful feature of the terminal is the **awk** programming
 language, which allows us to extract subsets of a data file, filter
 data according to some criteria or perform arithmetic operations on
 the data. awk manipulates a data file by per- forming operations on
 its **columns** - this is extremely useful for scientific data sets
 because typically the columns features or variables of interest.

 For example, we can use awk to produce a new results file that only
 contains SNP rsIDs (column 1), allele frequencies (column 4) and *P*
 -values (column 7) as follows:

      awk '{ print $1,$4,$7}' GIANT_Height.txt > GIANT_Height_3cols.txt

 We can also use a \"conditional statement\" in awk to extract all
 significant [SNPs]

 from the results file, using the following command:

     awk '{if($7 < 5e-8) { print } }' GIANT_Height.txt > Significant_SNPs.txt

Or the short form:

     awk '$7 < 5e-8{ print}' GIANT_Height.txt > Significant_SNPs.txt

 "if($7<5e-8)" and "$7 < 5e-8" tell awk to extract any rows
 with column 7 (the column containing *P* -value) with a value of
 smaller than 5e-8 and {print} means that we would like to print the
 entire row when this criterion is met.


# Introduction to R

**R** is a useful programming language that allows us to perform a
variety of statis- tical tests and data manipulation. It can also be
used to generate fantastic data visualisations. Here we will go
through some of the basics of **R** so that you can better understand
the practicals throughout the workshop.


## Basics

 If you are not using R Studio then you can type **R** in your terminal
 to run **R** in the terminal.

 ## Adding script to working dir

      cd ~/data/Day1a_Data/Day1a_Data
      wget https://raw.githubusercontent.com/WCSCourses/prs_2023/main/scripts/nagelkerke.R

 

## Working Directory

 When we start **R**, we will be working in a specific folder called
 the **working directory**. We can check the current/working
 directory we are in by typing:

       getwd()

 And we can change our working directory to the **Practical** folder by

       setwd("~/data/Day1a_Data/Day1a_Data")

## Libraries

 Most functionality of **R** is organised in \"packages\" or
 \"libraries\". To access these functions, we will have to install and
 \"load\" these packages. Most commonly used packages are installed
 together with the standard installation process. You can install a new
 library using the install.packages function.

 For example, to install *ggplot2*, run the command:

        install.packages("ggplot2")


 After installation, you can load the library by typing

        library(ggplot2)

 Alternatively, we can import functions (e.g. that we have written)
 from an R script file on our computer. For example, you can load the
 Nagelkerke *R2* function by typing

         source("nagelkerke.R")

 And you are now able to use the Nagelkerke *R2* function (we will use
 this function at the end of this worksheet).

## Variables in R

 You can assign a value or values to any variable you want using \<-.
 e.g

Assign a number to a

      a <- 1

Assign a vector containing a,b,c to v1

      v1 <- c("a", "b","c")

## Functions

 You can perform lots of operations in **R** using diﬀerent built-in R
 functions. Some examples are below:

Assign number of samples

        nsample <- 10000

Generate nsample random normal variable with mean = 0 and sd = 1

        normal <- rnorm(nsample, mean=0,sd=1)

        normal.2 <- rnorm(nsample, mean=0,sd=1)

We can examine the first few entries of the result using head

      head(normal)

And we can obtain the mean and sd using

      mean(normal)

      sd(normal)

We can also calculate the correlation between two variables
 using cor

      cor(normal, normal.2)

## Plotting

 While **R** contains many powerful plotting functions in its base
 packages,customisationn can be diﬃcult (e.g. changing the colour
 scales, arranging the axes). **ggplot2** is a powerful visualization package that provides extensive
 flexibility and customi- sation of plots. As an example, we can do the following

Load the package

      library(ggplot2)

Specify sample size

      nsample <-1000

 Generate random grouping using sample with replacement

     groups <- sample(c("a","b"), nsample, replace=T)

 Now generate the data

     dat <- data.frame(x=rnorm(nsample), y=rnorm(nsample), groups)

 Generate a scatter plot with diﬀerent coloring based on group

     ggplot(dat, aes(x=x,y=y,color=groups))+geom_point()


## Regression Models

 In statistical modelling, regression analyses are a set of statistical
 techniques for estimating the relationships among variables or
 features. We can perform regression analysis in **R**.

 Use the following code to perform linear regression on simulated
 variables "x" and "y":

 Simulate data

     nsample <- 10000

     x <- rnorm(nsample)

     y <- rnorm(nsample)

 Run linear regression

      lm(y~x)

 We can store the result into a variable

       reg <- lm(y~x)

 And get a detailed output using summary

       summary(lm(y~x))

 We can also extract the coeﬃcient of regression using

       reg$coeﬃcient

 And we can obtain the residuals by

      residual <- resid(reg)

 Examine the first few entries of residuals

       head(residual)

 We can also include covariates into the model

      covar <- rnorm(nsample)

      lm(y~x+covar)

 And can even perform interaction analysis

      lm(y~x+covar+x*covar)

 Alternatively, we can use the glm function to perform the regression:

      glm(y~x)

 For binary traits (case controls studies), logistic regression can be
 performed using

 Simulate samples

      nsample <- 10000

      x <- rnorm(nsample)

Simulate binary traits (must be coded with 0 and 1)

      y <- sample(c(0,1), size=nsample, replace=T)

 Perform logistic regression

      glm(y~x, family=binomial)

 Obtain the detailed output

      summary(glm(y~x, family=binomial))

 We will need the NagelkerkeR2 function

 to calculate the pseudo R2 for logistic model

       source("nagelkerke.R")

       reg <- glm(y~x, family=binomial)

 Calculate the Nagelkerke R2 using the NagelkerkeR2 function

         NagelkerkeR2(reg)
