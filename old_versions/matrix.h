/* 
 * File:   matrix.h
 * Author: fake
 *
 * Created on June 27, 2012, 3:55 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H

template <typename data_t>
class matrix {
    data_t *buffer;
    int nrows;
    int ncols;
    int total;
public:

    matrix(int rows, int cols, int fill=0) {
        total = rows*cols;
        nrows = rows;
        ncols = cols;
        if (nrows) {
            fprintf(stderr, "Requested allocation space of %d kbytes\n", (((int)sizeof (data_t)) * total) >> 10);
            buffer = new data_t[total];
            if(fill) 
                for (int i=0; i<total ;i++)
                        buffer[i]=fill;
        }
    }

    data_t & operator()(int row, int col) {
        if(row>=nrows || col>=ncols){
            fprintf(stderr,"Requesting element out of bounds at (%d, %d)!\n",row,col);
        }
        return buffer[ncols * row + col];
    }

    ~matrix() {
        if(nrows)
                delete buffer;
    }
};

#endif	/* MATRIX_H */

