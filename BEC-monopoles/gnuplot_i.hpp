////////////////////////////////////////////
//
// A C++ interface to gnuplot. 
//
// This is a direct translation from the C interface
// written by N. Devillard (which is available from
// http://ndevilla.free.fr/gnuplot/).
//
// As in the C interface this uses pipes and so wont
// run on a system that does'nt have POSIX pipe 
// support
//
// Rajarshi Guha
// <rajarshi@presidency.com>
//
// 07/03/03
//
// /////////////////////////////////////////

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

#include <stdarg.h>
#include <unistd.h>

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <stdexcept>

#define GP_MAX_TMP_FILES    64
#define GP_TMP_NAME_SIZE    512
#define GP_CMD_SIZE         1024
#define GP_TITLE_SIZE       80

using namespace std;

class GnuplotException : public runtime_error
{
    public:
        GnuplotException(const string &msg) : runtime_error(msg){}
};

class Gnuplot
{
    private:
        FILE            *gnucmd;
        string           pstyle;
        vector<string>   to_delete;
        int              nplots;
        bool             get_program_path(const string);
    public:
        Gnuplot();

        // set a style during construction
        Gnuplot(const string &);
        
        ~Gnuplot();

        // send a command to gnuplot
        void cmd(const string &, ...);

        // set line style
        void set_style(const string &);

        // set y and x axis labels
        void set_ylabel(const string &);
        void set_xlabel(const string &);

        // plot a single vector
        void plot_x(vector<double>, 
                const string & // title
                );

        // plot x,y pairs
        void plot_xy(vector<double>, vector<double>, 
                const string  & // title
                );

        // plot an equation supplied as a string
        void plot_equation(
                const string &, // equation 
                const string &  // title
                );

        // if multiple plots are present it will clear the plot area
        void reset_plot(void);
        
};

#endif
