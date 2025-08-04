#Make sure you have lftp installed
lftp -c "set ftp:list-options -a;
open https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/607/E-MTAB-10607/Files/;
mirror --verbose"
