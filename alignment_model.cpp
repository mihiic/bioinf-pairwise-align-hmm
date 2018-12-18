#include "alignment_model.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cfloat>

using namespace std;

AlignmentModel::AlignmentModel(unsigned long sample) {
    sampleSize = sample;
    transitions = {
            {0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0},   // t(0, 0), t(0, 1), t(0, 2), t(0, 3), t(0, 4)
            {0, 0.7,       0.1,       0.1,       0.1}, // t(1, 0), t(1, 1), t(1, 2), t(1, 3), t(1, 4)
            {0, 0.25,      0.65,      0,         0.1}, // t(2, 0), t(2, 1), t(2, 2), t(2, 3), t(2, 4)
            {0, 0.25,      0,         0.65,      0.1}, // t(3, 0), t(3, 1), t(3, 2), t(3, 3), t(3, 4)
            {0, 0,         0,         0,         0}    // t(4, 0), t(4, 1), t(4, 2), t(4, 3), t(4, 4)
    };

    emissions = {
            {0,    0.25,       0.25,       0.25,       0.25 /
                                                       1.0},  // {"-, -"}, {"-, A"}, {"-, T"}, {"-, C"}, {"-, G"},
            {0.25, 1.0 / 7.0,  1.0 / 28.0, 1.0 / 28.0, 1.00 /
                                                       28.0}, // {"A, -"}, {"A, A"}, {"A, T"}, {"A, C"}, {"A, G"},
            {0.25, 1.0 / 28.0, 1.0 / 7.0,  1.0 / 28.0, 1.0 /
                                                       28.0},  // {"T, -"}, {"T, A"}, {"T, T"}, {"T, C"}, {"T, G"},
            {0.25, 1.0 / 28.0, 1.0 / 28.0, 1.0 / 7.0,  1.0 /
                                                       28.0},  // {"C, -"}, {"C, A"}, {"C, T"}, {"C, C"}, {"C, G"},
            {0.25, 1.0 / 28.0, 1.0 / 28.0, 1.0 / 28.0, 1.0 /
                                                       7.0}    // {"G, -"}, {"G, A"}, {"G, T"}, {"G, C"}, {"G, G"},
    };

    currentLine = 0;
}

AlignmentModel::~AlignmentModel() {

}

double AlignmentModel::getEmissionProb(char i, char j) {
    return emissions[charToIndex(i)][charToIndex(j)];
}

int AlignmentModel::charToIndex(char x) {
    switch (x) {
        case '-':
            return 0;
        case 'A':
            return 1;
        case 'T':
            return 2;
        case 'C':
            return 3;
        case 'G':
            return 4;
        default:
            break;
    }
    return 0;
}

void AlignmentModel::loadTransition(char *filename) {
    ifstream file(filename);
    transitions.resize(0);
    transitions.resize(5, {});

    char value[256];
    char c;
    int col = 0;
    int row = 0;
    int index = 0;
    while (file.get(c)) {
        if (c == ',' || c == '\n') {
            value[index] = '\0';
            string val = string(value);
            transitions[row].push_back(stod(value, NULL));
            index = 0;
            col++;
            if (col == 5) {
                row++;
                col = 0;
            }
            continue;
        }
        if (c == ' ') {
            continue;
        }
        value[index] = c;
        index++;
    }

    file.close();
}

void AlignmentModel::loadEmission(char *filename) {
    ifstream file(filename);
    emissions.resize(0);
    emissions.resize(5, {});

    char value[256];
    char c;
    int col = 0;
    int row = 0;
    int index = 0;
    while (file.get(c)) {
        if (c == ',' || c == '\n') {
            value[index] = '\0';
            string val = string(value);
            emissions[row].push_back(stod(value, NULL));
            index = 0;
            col++;
            if (col == 5) {
                row++;
                col = 0;
            }
            continue;
        }
        if (c == ' ') {
            continue;
        }
        value[index] = c;
        index++;
    }

    file.close();
}

void AlignmentModel::run(char *filenameA, char *filenameB, char *outputFilename, char *transmissionFilename,
                         char *emissionFilename) { //
    loadTransition(transmissionFilename);
    loadEmission(emissionFilename);

    ifstream fileA(filenameA);
    ifstream fileB(filenameB);
    ofstream outA("predictions/a.fasta");
    ofstream outB("predictions/b.fasta");

    L1 = 0;
    L2 = 0;

    calculateDimensions(fileA, fileB);
    writePredictionHeaders(outA, outB);

    string lineA;
    string lineB;
    for (int i = 1; i <= totalLines; i++) {
        reserveMemory(maxL1, maxL2);
        lineA = readLine(fileA, (int) maxL1);
        lineB = readLine(fileB, (int) maxL2);
        double m = recursion(maxL1, maxL2, lineA, lineB);
        termination(maxL1, maxL2, m);
        traceback(maxL1, maxL2);
        convertPredictedAlignment(lineA, lineB);
        writePredictions(outA, outB);
    }

    reserveMemory(leftover1, leftover2);
    lineA = readLine(fileA, (int) leftover1);
    lineB = readLine(fileB, (int) leftover2);
    double m = recursion(leftover1, leftover2, lineA, lineB);
    termination(leftover1, leftover2, m);
    traceback(leftover1, leftover2);
    convertPredictedAlignment(lineA, lineB);
    writePredictions(outA, outB);

    fileA.close();
    fileB.close();
    outA.close();
    outB.close();

    concatePredictions(outputFilename);
}

void AlignmentModel::reserveMemory(unsigned long l1, unsigned long l2) {
    viterbi.resize(0);
    point.resize(0);

    viterbi.resize(l1 + 2, {});
    point.resize(l1 + 2, {});
    for (long i = 0; i < l1 + 2; i++) {
        viterbi[i].resize(l2 + 2, {});
        point[i].resize(l2 + 2, {});
        for (long j = 0; j < l2 + 2L; j++) {
            viterbi[i][j].resize(5, -DBL_MAX);
            point[i][j].resize(5, -1);
        }
    }

    viterbi[0][0][0] = 0;
}

void AlignmentModel::calculateDimensions(ifstream &fileA, ifstream &fileB) {
    char c;
    bool passedFirstLine = false;
    int i = 0;
    while (fileA.get(c)) {
        if (c == '\n' && !passedFirstLine) {
            passedFirstLine = true;
            header1[i] = '\0';
            continue;
        }
        if (passedFirstLine) {
            if (c == 'A' || c == 'G' || c == 'T' || c == 'C') {
                L1 += 1;
            }
        } else {
            header1[i] = c;
            i++;
        }
    }

    i = 0;
    passedFirstLine = false;
    while (fileB.get(c)) {
        if (c == '\n' && !passedFirstLine) {
            passedFirstLine = true;
            header2[i] = '\0';
            continue;
        }
        if (passedFirstLine) {
            if (c == 'A' || c == 'G' || c == 'T' || c == 'C') {
                L2 += 1;
            }
        } else {
            header2[i] = c;
            i++;
        }
    }

    cout << header1 << endl;
    cout << header2 << endl;

//    L2 -= 1;
//    L1 -= 1;

    double tl;
    if (L1 > L2) {
        maxL1 = sampleSize;
        totalLines = (L1 / maxL1);
        tl = (double) L1 / maxL1;
        maxL2 = (unsigned long) (L2 / tl);
    } else {
        maxL2 = sampleSize;
        totalLines = (L2 / maxL2);
        tl = (double) L2 / maxL2;
        maxL1 = (unsigned long) (L2 / tl);
    }

    leftover1 = L1 - maxL1 * totalLines;
    leftover2 = L2 - maxL2 * totalLines;

    fileA.clear();
    fileB.clear();
    fileA.seekg(0, ios::beg);
    fileB.seekg(0, ios::beg);

    // skip comment lines in both files
    while (fileA.get(c)) {
        if (c == '\n') {
            break;
        }
    }

    while (fileB.get(c)) {
        if (c == '\n') {
            break;
        }
    }
}

double AlignmentModel::recursion(long l1, long l2, string &line1, string &line2) {
    double max = -DBL_MAX;
    for (int i = 0; i <= l1; i++) {
        for (int j = 0; j <= l2; j++) {
            if (i == 0 && j == 0) {
                continue;
            }

            for (int k = 1; k <= 3; k++) {
                if (k == 1 && i > 0 && j > 0) {
                    for (int m = 0; m <= 4; m++) {
                        if (max < viterbi[i - 1][j - 1][m] + log(transitions[m][k])) {
                            max = viterbi[i - 1][j - 1][m] + log(transitions[m][k]);
                            viterbi[i][j][k] = log(getEmissionProb(
                                    line1.at((unsigned long) (i - 1)), line2.at((unsigned long) (j - 1))
                            )) + viterbi[i - 1][j - 1][m] + log(transitions[m][k]);

                            point[i][j][k] = m;
                        }
                    }
                    max = -DBL_MAX;
                } else if (k == 2 && i > 0) {
                    for (int m = 0; m <= 4; m++) {
                        if (max < viterbi[i - 1][j][m] + log(transitions[m][k])) {
                            max = viterbi[i - 1][j][m] + log(transitions[m][k]);
                            viterbi[i][j][k] = log(getEmissionProb(
                                    line1.at((unsigned long) i - 1), '-'
                            )) + viterbi[i - 1][j][m] + log(transitions[m][k]);
                            point[i][j][k] = m;
                        }
                    }
                    max = -DBL_MAX;
                } else if (k == 3 && j > 0) {
                    for (int m = 0; m <= 4; m++) {
                        if (max < viterbi[i][j - 1][m] + log(transitions[m][k])) {
                            max = viterbi[i][j - 1][m] + log(transitions[m][k]);
                            viterbi[i][j][k] = log(getEmissionProb(
                                    '-', line2.at((unsigned long) j - 1)
                            )) + viterbi[i][j - 1][m] + log(transitions[m][k]);
                            point[i][j][k] = m;
                        }
                    }
                    max = -DBL_MAX;
                }
            }
        }
    }
//    cout << max << endl;
    return max;
}

string AlignmentModel::readLine(ifstream &file, int lineSize) {
    char line[4096];
    int currentIndex = 0;
    char c;
    while (currentIndex < lineSize && file.get(c)) {
        if (c == '\n') {
            continue;
        }

        line[currentIndex] = c;
        currentIndex++;
    }
    line[currentIndex] = '\0';

    return string(line);
}

void AlignmentModel::termination(long l1, long l2, double max) {
    for (int m = 0; m <= 4; m++) {
        if (max < viterbi[l1][l2][m] + log(transitions[m][4])) {
            max = viterbi[l1][l2][m] + log(transitions[m][4]);
            viterbi[l1 + 1][l2 + 1][4] = viterbi[l1][l2][m] + log(transitions[m][4]);
            point[l1 + 1][l2 + 1][4] = m;
        }
    }
}

void AlignmentModel::traceback(long l1, long l2) {
    statePath.resize(0);
    statePath.push_back(4);
    int statePathIndex = 0;
    long l1Index = l1 + 1;
    long l2Index = l2 + 1;

    while ((l1Index > 0 || l2Index > 0) && statePath.at((unsigned long) statePathIndex) != 0) {
        long currState = statePath.at((unsigned long) statePathIndex);
        statePath.push_back(point[l1Index][l2Index][currState]);
        if (currState == 1 || currState == 4) {
            l1Index--;
            l2Index--;
        } else if (currState == 2) {
            l1Index--;
        } else if (currState == 3) {
            l2Index--;
        }
        statePathIndex++;
    }
}

void AlignmentModel::convertPredictedAlignment(string &lineA, string &lineB) {
    vector<char> statePathCharacters;
    for (long i = statePath.size() - 1; i >= 0; i--) {
        statePathCharacters.push_back((char) (statePath.at(i) + '0'));
    }

    predictedAlignmentA.resize(statePath.size() - 2);
    predictedAlignmentB.resize(statePath.size() - 2);

    long index1 = 0;
    long index2 = 0;

    for (int i = 1; i < statePathCharacters.size() - 1; i++) {
        if (statePathCharacters[i] == '1') {
            predictedAlignmentA[i - 1] = lineA[index1];
            predictedAlignmentB[i - 1] = lineB[index2];
            index1++;
            index2++;
        } else if (statePathCharacters[i] == '2') {
            predictedAlignmentA[i - 1] = lineA[index1];
            predictedAlignmentB[i - 1] = '-';
            index1++;
        } else if (statePathCharacters[i] == '3') {
            predictedAlignmentA[i - 1] = '-';
            predictedAlignmentB[i - 1] = lineB[index2];
            index2++;
        }
    }
}

void AlignmentModel::writePredictionHeaders(ofstream &fileA, ofstream &fileB) {
    for (char c : header1) {
        if (c == '\0') {
            break;
        }
        fileA << c;
    }
    for (char c : header2) {
        if (c == '\0') {
            break;
        }
        fileB << c;
    }

    fileA << endl;
    fileB << endl;
    fileA.flush();
    fileB.flush();
}

void AlignmentModel::writePredictions(ofstream &fileA, ofstream &fileB) {
    int temp = currentLine;
    for (char i : predictedAlignmentA) {
        fileA << i;
        temp++;
        if (temp > 60) {
            temp = 0;
            fileA << endl;
        }
    }

    temp = currentLine;
    for (char i : predictedAlignmentB) {
        fileB << i;
        temp++;
        if (temp > 60) {
            temp = 0;
            fileB << endl;
        }
    }

    currentLine = temp;
    fileA.flush();
    fileB.flush();
}

void AlignmentModel::concatePredictions(char *filename) {
    ifstream file1("predictions/a.fasta");
    ifstream file2("predictions/b.fasta");
    ofstream combined_file(filename);
    combined_file << file1.rdbuf() << endl << file2.rdbuf();

    combined_file.close();
    file1.close();
    file2.close();

    remove("predictions/a.fasta");
    remove("predictions/b.fasta");
}

