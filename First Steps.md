# Task to do
* Make a github repo (done)
* Copy the transcriptome data to Xue's info account
* Practce using Screen
* Install Trinity
* Run assemblies for XL and XT
* Download Xl and XT genomes
* Make BLAST databases from these genomes.


# Basics for Unix

If you want to cancel a command type `Command .`

If you want to view the contents of a directory type `ls`.

If you want to change a directory type `cd `\<directory path\>.

If you want to know what directory you are in type `pwd`.

If you want to make a directory, type `mkdir` \<directory name\>.

If you want to remove a file type `rm -f `\<filename\>.

To download the trinity package we typed `wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.2.0.tar.gz`.  This link would have to be updated when an new version becomes available.

To uncompress it, type `tar -zxvf v2.2.0.tar.gz`

Instructions for install are her `https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing%20Trinity`.

# Basics for using info.

To connect to info, use this command:

`ssh xue@info.mcmaster.ca`

This gets you into the `head` server, which is called infoserv.  Please do not run jobs on the head.  Instead you can connect to the other processors like this: `rsh info115` (or `rsh info113`, `rsh info114`,...).

You can check how busy a system is using the unix command `top`.  Exit by typing `q`.

# Using Unix Screen

Screen is very useful for running background jobs.

To launch a screen session type `screen -S `\<session name\>.

To exit a screen type `control a, control d`.  

To list screen type `screen -ls`.

To reconnect with a screen type `screen -r `\<session name\>.

To kill a screen type `screen -X -S ` \<session name\> ` kill`.

To exit the server type `logout`.
