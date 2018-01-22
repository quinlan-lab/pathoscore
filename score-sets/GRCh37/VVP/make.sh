if [ ! -s VVP_scores_allChr.txt.gz ]; then
    wget http://weatherby.genetics.utah.edu/VVP_Scores/VVP_scores_allChr.txt.gz
    wget http://weatherby.genetics.utah.edu/VVP_Scores/VVP_scores_allChr.txt.gz.tbi
fi
