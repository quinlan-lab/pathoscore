#use 11th column (MTR)
if [ ! -s mtrflatfile_1.0.txt.gz ]; then
    wget http://mtr-viewer.mdhs.unimelb.edu.au:8079/mtrflatfile_1.0.txt.gz
    wget http://mtr-viewer.mdhs.unimelb.edu.au:8079/mtrflatfile_1.0.txt.gz.tbi
fi
