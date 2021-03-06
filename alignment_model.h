#ifndef BIONF_ALIGNMENT_MODEL_H
#define BIONF_ALIGNMENT_MODEL_H

#include <vector>
#include <fstream>

using namespace std;

class AlignmentModel {
public:
    AlignmentModel(unsigned long);
    ~AlignmentModel();

    vector<vector<double>> transitions;
    vector<vector<double>> emissions;
    unsigned long sampleSize;

    double getEmissionProb(char i, char j);
    int charToIndex(char x);

    void run(char* filenameA, char* filenameB, char* outputFilename, char* transmissionFilename, char* emissionFilename); //
    void loadTransition(char* filename);
    void loadEmission(char* filename);
    void reserveMemory(unsigned long, unsigned long);
    void calculateDimensions(ifstream &fileA, ifstream &fileB);
    double recursion(long, long, string&, string&);
    string readLine(ifstream&, int);

    void termination(long, long, double);
    void traceback(long, long);
    void convertPredictedAlignment(string&, string&);

    void writePredictionHeaders(ofstream&, ofstream&);
    void writePredictions(ofstream&, ofstream&);

    void concatePredictions(char* filename);

    unsigned long L1;
    unsigned long L2;
    unsigned long maxL1;
    unsigned long maxL2;
    unsigned long leftover1;
    unsigned long leftover2;
    unsigned long totalLines;

    int currentLine;

    char header1[256];
    char header2[256];

    vector<vector<vector<double>>> viterbi;
    vector<vector<vector<long>>> point;

    vector<long> statePath;
    vector<char> predictedAlignmentA;
    vector<char> predictedAlignmentB;

    static const int BLANK = 0;
    static const int A = 1;
    static const int B = 2;
    static const int C = 3;
    static const int D = 4;
};


#endif