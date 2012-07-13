/* 
 * File:   matrix.h
 * Author: fake
 *
 * Created on June 27, 2012, 3:55 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H

#include <vector>

template <typename data_t>
class matrix {
    std::vector<data_t> buffer;
    int nrows;
    int ncols;
    int total;
public:

    matrix(int rows, int cols, data_t fill = 0) {
        total = rows*cols;
        nrows = rows;
        ncols = cols;
        if (nrows) {
#ifdef DEBUG
            fprintf(stderr, "Requested allocation space of %d kbytes\n", (((int) sizeof (data_t)) * total) >> 10);
#endif
            buffer.reserve(total);
            buffer.assign(total, fill);
        }
    }

    data_t & operator()(int row, int col) {
        if (row >= nrows || col >= ncols) {
            //fprintf(stderr,"Requesting element out of bounds at (%d, %d)!\n",row,col);
        }
        return buffer[ncols * row + col];
    }
};

template <int cols>
class memory {
    std::vector<int> buffer;
    int nrows;
    int total;
public:

    memory(int rows = 0, int fill = 0) {
        total = rows*cols;
        nrows = rows;
        if (nrows) {
#ifdef DEBUG
            fprintf(stderr, "Requested allocation space of %d kbytes\n", (((int) sizeof (int)) * total) >> 10);
#endif
            buffer.reserve(total);
            buffer.assign(total, fill);
        }
    }

    int & operator()(int row, int col) {
        return buffer[cols * row + col];
    }

    int & operator[](int row) {
        return buffer[row];
    }

    int & operator()(int row) {
        return buffer[cols * row];
    }
};

#endif	/* MATRIX_H */

