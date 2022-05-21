// algorithm from https://github.com/orhancelik1/Needleman-Wunsh-Algorithm

#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include "kseq.h"
#include <zlib.h>  
#include <iostream>
#include <vector>

using namespace std;

int match = 5, gap = -5, missmatch = -3; // scoring line

KSEQ_INIT(gzFile, gzread)  



int Max(int num1, int num2, int num3)
{
    int max;
    if (num1 > num2) {
        if (num1 > num3) {
            max = num1;
        }
        else {
            max = num3;
        }
    } else {
        if (num2 > num3) {
            max = num2;
        } else {
            max = num3;
        }
    }
    return max;
}

int cal_indentity(const char *seq1, const char *seq2)
{
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    int maxlenori = Max(len1, len2, 0);
    int i, j, num1, num2, num3, c1 = 0, c2 = 0;
    char Tbseq1[2*maxlenori]; // traceback array for seq1
    char Tbseq2[2*maxlenori]; // traceback array for seq2
    int lenseq1 = strlen(seq1);
    int lenseq2 = strlen(seq2);
    int scoreMatrix[lenseq1][lenseq2];

    // building matrix and putting 'alpha' values
    // printf("\t");
    // for(j=0;j<lenseq2;j++){
    //	printf("%c \t",seq2[j]);
    //}
    // printf("\n");

    for (i = 0; i < lenseq1; i++) {
        // printf("%c \t",seq1[i]);
        for (j = 0; j < lenseq2; j++) {
            if (i == 0 || j == 0) { // gap penalty values
                if (i == 0) {
                    scoreMatrix[0][j] = j * gap;
                    // printf("%d\t",scoreMatrix[0][j]);
                }
                else { // j==0
                    scoreMatrix[i][0] = i * gap;
                    // printf("%d\t",scoreMatrix[i][0]);
                }
            }
            else {
                if (seq1[i] == seq2[j]) {
                    num1 = scoreMatrix[i - 1][j - 1] + match;
                    num2 = scoreMatrix[i - 1][j] + gap;
                    num3 = scoreMatrix[i][j - 1] + gap;
                    // printf("%d\t",scoreMatrix[i][j]=Max(num1,num2,num3));
                    scoreMatrix[i][j] = Max(num1, num2, num3);
                } else {
                    num1 = scoreMatrix[i - 1][j - 1] + missmatch;
                    num2 = scoreMatrix[i - 1][j] + gap;
                    num3 = scoreMatrix[i][j - 1] + gap;
                    // printf("%d\t",scoreMatrix[i][j]=Max(num1,num2,num3));
                    scoreMatrix[i][j] = Max(num1, num2, num3);
                }
            }
        }
        // printf("\n\n");
    }
    // traceback part

    j = lenseq2 - 1, i = lenseq1 - 1;
    while (i > 0)
    {
        if (seq1[i] == seq2[j])
        { // match condition
            Tbseq1[c1] = seq1[i]; // Tbseq traceback sequence
            Tbseq2[c2] = seq2[j];
            c1++, c2++; // Tbseq indises
            i--, j--;	// diagonal move for match
        }
        else
        {
            if (scoreMatrix[i - 1][j - 1] + missmatch == scoreMatrix[i][j])
            { // missmatch condition

                Tbseq1[c1] = seq1[i];
                Tbseq2[c2] = seq2[j];
                c1++, c2++;
                i--, j--; // diagonal move for missmatch
            }

            else if (scoreMatrix[i][j - 1] + gap == scoreMatrix[i][j])
            { // gap condition

                Tbseq1[c1] = '-';
                Tbseq2[c2] = seq2[j];
                c1++, c2++;
                j--; // move left
            }

            else if (scoreMatrix[i - 1][j] + gap == scoreMatrix[i][j])
            { // gap condition

                Tbseq1[c1] = seq1[i];
                Tbseq2[c2] = '-';
                c1++, c2++;
                i--; // move top
            }
        }
    }

    //printf("Optimal Global Alignments:\n");
    //printf("sequence 1:");
    for (i = lenseq1 - 1; i >= 0; i--) {
        //printf("%c", Tbseq1[i]);
    }
    //printf("\n");
    //printf("sequence 2:");
    for (i = lenseq1 - 1; i >= 0; i--) {
        //printf("%c", Tbseq2[i]);
    }
    //printf("\n");
    int same = 0;
    int all = 0;
    for (i = lenseq1 - 1; i >= 0; i--) {
        if (Tbseq1[i] != '-') {
            all++;
            // check is same
            if (Tbseq1[i] == Tbseq2[i]) {
            same++;
            }
        }
    }
    int identidy = 100 * same / all;
    return(same);
}





vector<string*> readseqs(gzFile *fp) {
    kseq_t *seq;
    seq = kseq_init(*fp); // STEP 3: initialize seq  
    int l;  
    vector<string*> reads;
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence  
        //printf("name: %s\n", seq->name.s);  
        //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);  
        //printf("seq: %s\n", seq->seq.s);  
        //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        string *strPtr=new string;
        strPtr->assign(seq->seq.s, seq->seq.l);
        reads.push_back(strPtr);
    }
    kseq_destroy(seq); // STEP 5: destroy seq  
    return(reads);
}

int main(int argc, char *argv[])  
{
    gzFile fp;
    // open input file from stdin if no file name given
    if (argc == 1) {
        fp = gzdopen(fileno(stdin), "r");  
    } else {  
        fp = gzopen(argv[1], "r");   // STEP 2: open the file handler  
    }
    vector<string*> seqs = readseqs(&fp);
    int nseqs = seqs.size();

    // loop seqs pairs
    for(int i1=0; i1<nseqs; i1++)
    {
        for(int i2=i1+1; i2<nseqs; i2++)
        {
            int len1 = seqs[i1]->length();
            int len2 = seqs[i2]->length();
            string seq1n = "$";
            seq1n = seq1n.append( *seqs[i1] );
            string seq2n = "$";
            seq2n = seq2n.append( *seqs[i1] );
            //cout << seq1n << endl;
            //cout << seq2n << endl;
            float meanlen = (len1 + len2) / 2;
            int same = cal_indentity(seqs[i1]->c_str(), seqs[i2]->c_str());
            float identity = 100 * same / meanlen;
            printf("%d\t%d\t%f\t%d\t%d\t%d\n", i1, i2, identity,  len1, len2, same );
        }
    }
    return 0;  
}




