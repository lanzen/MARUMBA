import os

for file in os.listdir("/share/raw/MiSeq/2021/MARUMBA"):
    if file.endswith("R1_001.fastq.gz"):        
        stub = file[file.find("-")+1:file.find("_")]
        os.mkdir(stub)
#        os.chdir(stub)
        os.symlink("/share/raw/MiSeq/2021/MARUMBA/"+file, stub+"/"+file)
        r2 = file.replace("R1_001.fastq.gz","R2_001.fastq.gz")
        os.symlink("/share/raw/MiSeq/2021/MARUMBA/"+r2, stub+"/"+r2)
#        os.chdir("..")
        
