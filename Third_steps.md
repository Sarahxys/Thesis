# Testing out Xue's new script!

First, we need to make sure that the script is executable like this:

```
chmod 755 BlastnOutputAnalysis.pl
```

To execute the script locally (not on server):
```
./BlastnOutputAnalysis.pl
```
To put perl script on the server so that I can run it on the server:
```
scp BlastnOutputAnanlysis_Oct24_v3_clean.pl  xue@info.mcmaster.ca:/home/xue/transcriptome_data/blastn_BJE3909_output
```

#To empty space in head 1 server, move it and zip it 
To know how large the file is:
 ```
 du -sh *
 ```
To move it to head 4:
```
mv trinityrnaseq-2.2.0/ /4/xue/
```
To create a symbolic link:
```
ln -s /4/xue/trinityrnaseq-2.2.0/ hereistrin
```
To move it back to head 1 from head 4:
```
mv /4/xue/trinityrnaseq-2.2.0/ .
```
To zip it:
```
gzip filename
```



