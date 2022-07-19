
git config --global user.email wangmeijiao2008@hotmail.com 
git config --global user.name wangmeijiao
git config --global init.defaultBranch main



##init a new repository
echo "# Human_placenta_multi-omics_snRNA-seq_snATAC-seq" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/wangmeijiao/Human_placenta_multi-omics_snRNA-seq_snATAC-seq.git
git push -u origin main


#token
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



##add files (add commit push)

git status #check 

git add . #add current dir
git commit -a -m "commit" #message
git push -u origin main




