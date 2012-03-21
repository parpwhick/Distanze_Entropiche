#include <string>
using namespace std;

void translate(string type, char *testo,int len) {
    char translation_table[128];

    for(int i=0; i<128;i++)
        translation_table[i]=i;
    
    translation_table['X']='.';
 
    
    if (type == "2") { // PH
        translation_table['A'] = 'P';
        translation_table['G'] = 'P';
        translation_table['T'] = 'P';
        translation_table['S'] = 'P';
        translation_table['N'] = 'P';
        translation_table['Q'] = 'P';
        translation_table['D'] = 'P';
        translation_table['E'] = 'P';
        translation_table['H'] = 'P';
        translation_table['R'] = 'P';
        translation_table['K'] = 'P';
        translation_table['P'] = 'P';
        translation_table['C'] = 'H';
        translation_table['M'] = 'H';
        translation_table['F'] = 'H';
        translation_table['I'] = 'H';
        translation_table['L'] = 'H';
        translation_table['V'] = 'H';
        translation_table['W'] = 'H';
        translation_table['Y'] = 'H';

    }
    if (type == "5") { // ARCTD
        translation_table['I'] = 'A';
        translation_table['V'] = 'A';
        translation_table['L'] = 'A';
        translation_table['F'] = 'R';
        translation_table['Y'] = 'R';
        translation_table['W'] = 'R';
        translation_table['H'] = 'R';
        translation_table['K'] = 'C';
        translation_table['R'] = 'C';
        translation_table['D'] = 'C';
        translation_table['E'] = 'C';
        translation_table['G'] = 'T';
        translation_table['A'] = 'T';
        translation_table['C'] = 'T';
        translation_table['S'] = 'T';
        translation_table['T'] = 'D';
        translation_table['M'] = 'D';
        translation_table['Q'] = 'D';
        translation_table['N'] = 'D';
        translation_table['P'] = 'D';

    }
    if (type == "6") { // ARPNTD
        translation_table['I'] = 'A';
        translation_table['V'] = 'A';
        translation_table['L'] = 'A';
        translation_table['F'] = 'R';
        translation_table['Y'] = 'R';
        translation_table['W'] = 'R';
        translation_table['H'] = 'R';
        translation_table['K'] = 'P';
        translation_table['R'] = 'P';
        translation_table['D'] = 'N';
        translation_table['E'] = 'N';
        translation_table['G'] = 'T';
        translation_table['A'] = 'T';
        translation_table['C'] = 'T';
        translation_table['S'] = 'T';
        translation_table['T'] = 'D';
        translation_table['M'] = 'D';
        translation_table['Q'] = 'D';
        translation_table['N'] = 'D';
        translation_table['P'] = 'D';

    }
    if (type == "3IMG") { // PNH
        translation_table['D'] = 'P';
        translation_table['N'] = 'P';
        translation_table['E'] = 'P';
        translation_table['Q'] = 'P';
        translation_table['K'] = 'P';
        translation_table['R'] = 'P';
        translation_table['G'] = 'N';
        translation_table['T'] = 'N';
        translation_table['S'] = 'N';
        translation_table['Y'] = 'N';
        translation_table['P'] = 'N';
        translation_table['M'] = 'N';
        translation_table['I'] = 'H';
        translation_table['V'] = 'H';
        translation_table['L'] = 'H';
        translation_table['F'] = 'H';
        translation_table['C'] = 'H';
        translation_table['M'] = 'H';
        translation_table['A'] = 'H';
        translation_table['W'] = 'H';
    }
    if (type == "5IMG") { // GCEMF  (IMGT amino acid volume)
        translation_table['G'] = 'G';
        translation_table['A'] = 'G';
        translation_table['S'] = 'G';
        translation_table['C'] = 'C';
        translation_table['D'] = 'C';
        translation_table['P'] = 'C';
        translation_table['N'] = 'C';
        translation_table['T'] = 'C';
        translation_table['E'] = 'E';
        translation_table['V'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['H'] = 'E';
        translation_table['M'] = 'M';
        translation_table['I'] = 'M';
        translation_table['L'] = 'M';
        translation_table['K'] = 'M';
        translation_table['R'] = 'M';
        translation_table['F'] = 'F';
        translation_table['Y'] = 'F';
        translation_table['W'] = 'F';

    }
    if (type == "11IMG") { // AFCGSWYPDNH
        translation_table['A'] = 'A';
        translation_table['V'] = 'A';
        translation_table['I'] = 'A';
        translation_table['L'] = 'A';
        translation_table['F'] = 'F';
        translation_table['C'] = 'C';
        translation_table['M'] = 'C';
        translation_table['G'] = 'G';
        translation_table['S'] = 'S';
        translation_table['T'] = 'S';
        translation_table['W'] = 'W';
        translation_table['Y'] = 'Y';
        translation_table['P'] = 'P';
        translation_table['D'] = 'D';
        translation_table['E'] = 'D';
        translation_table['N'] = 'N';
        translation_table['Q'] = 'N';
        translation_table['H'] = 'H';
        translation_table['K'] = 'H';
        translation_table['R'] = 'H';

    }
    if (type == "Murphy15") { // LCAGSTPFWEDNQKH
        translation_table['L'] = 'L';
        translation_table['V'] = 'L';
        translation_table['I'] = 'L';
        translation_table['M'] = 'L';
        translation_table['C'] = 'C';
        translation_table['A'] = 'A';
        translation_table['G'] = 'G';
        translation_table['S'] = 'S';
        translation_table['T'] = 'T';
        translation_table['P'] = 'P';
        translation_table['F'] = 'F';
        translation_table['Y'] = 'F';
        translation_table['W'] = 'W';
        translation_table['E'] = 'E';
        translation_table['D'] = 'D';
        translation_table['N'] = 'N';
        translation_table['Q'] = 'Q';
        translation_table['K'] = 'K';
        translation_table['R'] = 'K';
        translation_table['H'] = 'H';

    }
    if (type == "Murphy10") { // LCAGSPFEKH
        translation_table['L'] = 'L';
        translation_table['V'] = 'L';
        translation_table['I'] = 'L';
        translation_table['M'] = 'L';
        translation_table['C'] = 'C';
        translation_table['A'] = 'A';
        translation_table['G'] = 'G';
        translation_table['S'] = 'S';
        translation_table['T'] = 'S';
        translation_table['P'] = 'P';
        translation_table['F'] = 'F';
        translation_table['Y'] = 'F';
        translation_table['W'] = 'F';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['N'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['K'] = 'K';
        translation_table['R'] = 'K';
        translation_table['H'] = 'H';

    }
    if (type == "Murphy8") { // LASPFEKH
        translation_table['L'] = 'L';
        translation_table['V'] = 'L';
        translation_table['I'] = 'L';
        translation_table['M'] = 'L';
        translation_table['C'] = 'L';
        translation_table['A'] = 'A';
        translation_table['G'] = 'A';
        translation_table['S'] = 'S';
        translation_table['T'] = 'S';
        translation_table['P'] = 'P';
        translation_table['F'] = 'F';
        translation_table['Y'] = 'F';
        translation_table['W'] = 'F';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['N'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['K'] = 'K';
        translation_table['R'] = 'K';
        translation_table['H'] = 'H';

    }
    if (type == "Murphy4") { // LAFE
        translation_table['L'] = 'L';
        translation_table['V'] = 'L';
        translation_table['I'] = 'L';
        translation_table['M'] = 'L';
        translation_table['C'] = 'L';
        translation_table['A'] = 'A';
        translation_table['G'] = 'A';
        translation_table['S'] = 'A';
        translation_table['T'] = 'A';
        translation_table['P'] = 'A';
        translation_table['F'] = 'F';
        translation_table['Y'] = 'F';
        translation_table['W'] = 'F';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['N'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['K'] = 'E';
        translation_table['R'] = 'E';
        translation_table['H'] = 'E';

    }
    if (type == "Murphy2") { //PE
        translation_table['L'] = 'P';
        translation_table['V'] = 'P';
        translation_table['I'] = 'P';
        translation_table['M'] = 'P';
        translation_table['C'] = 'P';
        translation_table['A'] = 'P';
        translation_table['G'] = 'P';
        translation_table['S'] = 'P';
        translation_table['T'] = 'P';
        translation_table['P'] = 'P';
        translation_table['F'] = 'P';
        translation_table['Y'] = 'P';
        translation_table['W'] = 'P';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['N'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['K'] = 'E';
        translation_table['R'] = 'E';
        translation_table['H'] = 'E';

    }
    if (type == "Wang5") { // IAGEK
        translation_table['C'] = 'I';
        translation_table['M'] = 'I';
        translation_table['F'] = 'I';
        translation_table['I'] = 'I';
        translation_table['L'] = 'I';
        translation_table['V'] = 'I';
        translation_table['W'] = 'I';
        translation_table['Y'] = 'I';
        translation_table['A'] = 'A';
        translation_table['T'] = 'A';
        translation_table['H'] = 'A';
        translation_table['G'] = 'G';
        translation_table['P'] = 'G';
        translation_table['D'] = 'E';
        translation_table['E'] = 'E';
        translation_table['S'] = 'K';
        translation_table['N'] = 'K';
        translation_table['Q'] = 'K';
        translation_table['R'] = 'K';
        translation_table['K'] = 'K';

    }
    if (type == "Wang5v") { // ILAEK
        translation_table['C'] = 'I';
        translation_table['M'] = 'I';
        translation_table['F'] = 'I';
        translation_table['I'] = 'I';
        translation_table['L'] = 'L';
        translation_table['V'] = 'L';
        translation_table['W'] = 'L';
        translation_table['Y'] = 'L';
        translation_table['A'] = 'A';
        translation_table['T'] = 'A';
        translation_table['G'] = 'A';
        translation_table['S'] = 'A';
        translation_table['N'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['D'] = 'E';
        translation_table['E'] = 'E';
        translation_table['H'] = 'K';
        translation_table['P'] = 'K';
        translation_table['R'] = 'K';
        translation_table['K'] = 'K';

    }
    if (type == "Wang3") { // IAE
        translation_table['C'] = 'I';
        translation_table['M'] = 'I';
        translation_table['F'] = 'I';
        translation_table['I'] = 'I';
        translation_table['L'] = 'I';
        translation_table['V'] = 'I';
        translation_table['W'] = 'I';
        translation_table['Y'] = 'I';
        translation_table['A'] = 'A';
        translation_table['T'] = 'A';
        translation_table['H'] = 'A';
        translation_table['G'] = 'A';
        translation_table['P'] = 'A';
        translation_table['R'] = 'A';
        translation_table['D'] = 'E';
        translation_table['E'] = 'E';
        translation_table['S'] = 'E';
        translation_table['N'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['K'] = 'E';

    }
    if (type == "Wang2") { // IA
        translation_table['C'] = 'I';
        translation_table['M'] = 'I';
        translation_table['F'] = 'I';
        translation_table['I'] = 'I';
        translation_table['L'] = 'I';
        translation_table['V'] = 'I';
        translation_table['W'] = 'I';
        translation_table['Y'] = 'I';
        translation_table['A'] = 'A';
        translation_table['T'] = 'A';
        translation_table['H'] = 'A';
        translation_table['G'] = 'A';
        translation_table['P'] = 'A';
        translation_table['R'] = 'A';
        translation_table['D'] = 'A';
        translation_table['E'] = 'A';
        translation_table['S'] = 'A';
        translation_table['N'] = 'A';
        translation_table['Q'] = 'A';
        translation_table['K'] = 'A';

    }
    if (type == "Li10") { // CYLVGPSNEK
        translation_table['C'] = 'C';
        translation_table['F'] = 'Y';
        translation_table['Y'] = 'Y';
        translation_table['W'] = 'Y';
        translation_table['M'] = 'L';
        translation_table['L'] = 'L';
        translation_table['I'] = 'V';
        translation_table['V'] = 'V';
        translation_table['G'] = 'G';
        translation_table['P'] = 'P';
        translation_table['A'] = 'S';
        translation_table['T'] = 'S';
        translation_table['S'] = 'S';
        translation_table['N'] = 'N';
        translation_table['H'] = 'N';
        translation_table['Q'] = 'E';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['R'] = 'K';
        translation_table['K'] = 'K';
    }
    if (type == "Li5") { // YIGSE
        translation_table['C'] = 'Y';
        translation_table['F'] = 'Y';
        translation_table['Y'] = 'Y';
        translation_table['W'] = 'Y';
        translation_table['M'] = 'I';
        translation_table['L'] = 'I';
        translation_table['I'] = 'I';
        translation_table['V'] = 'I';
        translation_table['G'] = 'G';
        translation_table['P'] = 'S';
        translation_table['A'] = 'S';
        translation_table['T'] = 'S';
        translation_table['S'] = 'S';
        translation_table['N'] = 'E';
        translation_table['H'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['R'] = 'E';
        translation_table['K'] = 'E';
    }
    if (type == "Li4") { // YISE
        translation_table['C'] = 'Y';
        translation_table['F'] = 'Y';
        translation_table['Y'] = 'Y';
        translation_table['W'] = 'Y';
        translation_table['M'] = 'I';
        translation_table['L'] = 'I';
        translation_table['I'] = 'I';
        translation_table['V'] = 'I';
        translation_table['G'] = 'S';
        translation_table['P'] = 'S';
        translation_table['A'] = 'S';
        translation_table['T'] = 'S';
        translation_table['S'] = 'S';
        translation_table['N'] = 'E';
        translation_table['H'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['R'] = 'E';
        translation_table['K'] = 'E';
    }
    if (type == "Li3") { // ISE
        translation_table['C'] = 'I';
        translation_table['F'] = 'I';
        translation_table['Y'] = 'I';
        translation_table['W'] = 'I';
        translation_table['M'] = 'I';
        translation_table['L'] = 'I';
        translation_table['I'] = 'I';
        translation_table['V'] = 'I';
        translation_table['G'] = 'S';
        translation_table['P'] = 'S';
        translation_table['A'] = 'S';
        translation_table['T'] = 'S';
        translation_table['S'] = 'S';
        translation_table['N'] = 'E';
        translation_table['H'] = 'E';
        translation_table['Q'] = 'E';
        translation_table['E'] = 'E';
        translation_table['D'] = 'E';
        translation_table['R'] = 'E';
        translation_table['K'] = 'E';
    }

    for(int i=0;i<len;i++)
        if(testo[i]>='a' && testo[i]<='z')
            testo[i]+='A'-'a';
    for(int i=0;i<len;i++)
        testo[i]=translation_table[(int)testo[i]];
}
