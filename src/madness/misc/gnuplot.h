#ifndef MADNESS_GNUPLOT_H__INCUDED
#define MADNESS_GNUPLOT_H__INCUDED

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <string>
#include <vector>

namespace madness {
    class Gnuplot {
        FILE *f;     // pipe connection to gnuplot process
        pid_t pid;   // pid of gnuplot process
        FILE *ftee;  // filestream for data tee'd from gnuplot process
        
        // base case for unpacking datablock value
        template <typename T>
        void dbvalue(size_t i, const T& t) {
            char buf[256];
            snprintf(buf,sizeof(buf),"%16.8e",double(t[i]));
            (*this)(buf); // newline
        }
        
        // recursion case for unpacking datablock value
        template <typename T, typename... Ts>
        void dbvalue(size_t i, const T& t, Ts... values) {
            char buf[256];
            snprintf(buf,sizeof(buf),"%16.8e ",double(t[i])); // space
            (*this)(buf,false);
            dbvalue(i,values...);
        }
        
        // base case for unpacking plot value
        template <int n, typename T>
        void doplot(const char* name, const T& value0) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s using 1:%d", name, n);
            (*this)(buf);
        }
        
        // recursion case for unpacking plot value
        template <int n, typename T, typename... Ts>
        void doplot(const char* name, const T& value0, Ts... values) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s using 1:%d, ", name, n);
            (*this)(buf,false);
            
            doplot<n+1,Ts...>(name, values...);
        }
        
    public:
        Gnuplot& operator=(Gnuplot&&) = delete;
        Gnuplot& operator=(const Gnuplot&) = delete;
        Gnuplot(const Gnuplot&) = delete;
        Gnuplot(const std::string& cmd = "", const std::string& teefile = "") : f(0), ftee(0) {
            int p[2];
            if (pipe (p)) {
                throw "Pipe failed.";
            }        
            pid = fork ();
            if (pid ==  0) { // Child process.
                close(p[1]);
                dup2(p[0],STDIN_FILENO);
                close(p[0]);
                if (execlp ("gnuplot", "gnuplot", "-persist", NULL) == -1) {
                    //if (execlp ("cat", "cat", "-", NULL) == -1) {
                    fprintf(stderr,"Gnuplot: execlp failed for gnuplot ... plotting disabled\n");
                    exit(1);
                }
            }
            else if (pid < (pid_t) 0)  { // Failure
                throw "Fork failed.";
            }
            else { // Parent process
                close (p[0]);
                f = fdopen (p[1], "w");
            }
            if (teefile.size() > 0) {
                ftee = fopen(teefile.c_str(),"w");
                if (!ftee) {
                    fprintf(stderr,"Gnuplot: fopen failed for tee file %s ... tee of plotting disabled\n",teefile.c_str());
                }
            }
            
            if (cmd.size() > 0) (*this)(cmd);
        }
        
        // outputs string to gnuplot process
        void operator()(const char* cmd, bool EOL=true) {
            
            if (f) {
                if (!fprintf(f,"%s",cmd)) {
                    fprintf(stderr,"Gnuplot: failed writing to gnuplot pipe ... plotting disabled\n");
                    fclose(f);
                    f = NULL;
                }
            }
            if (ftee) fprintf(ftee,"%s",cmd);

            const int n = strlen(cmd);
            if (EOL && ((n==0) || (cmd[n-1] != '\n') ) ) {
                if (f) {
                    if (!fprintf(f,"\n")) {
                        fprintf(stderr,"Gnuplot: failed writing newline to gnuplot pipe ... plotting disabled\n");
                        fclose(f);
                        f = NULL;
                    }
                }
                if (ftee) fprintf(ftee,"\n");
            }
            if (f) fflush(f);
        }
        
        // outputs string to gnuplot process
        void operator()(const std::string& cmd, bool EOL=true) {
            (*this)(cmd.c_str(), EOL);
        }
        
        // Define a gnuplot data block with given name assuming 1-d indexing via [] with explicit size
        template <typename T, typename... Ts>
        void db(const std::string& name, size_t size, const T& x, Ts... values) {
            (*this)("$",false);
            (*this)(name,false);
            (*this)(" << EOD");
            for (size_t i = 0; i<size; ++i) {
                dbvalue(i,x,values...);
            }
            (*this)("EOD");
        }
        
        // Define a gnuplot data block with given name assuming 1-d indexing via [] with size from x (a 1-d container that supports size())
        template <typename T, typename... Ts>
        void db(const std::string& name, const T& x, Ts... values) {
            db(name,(size_t) x.size(),x,values...); // have to force x.size() to be size_t since Tensor<T>::size() returns long
        }
        
        // Plots data in 2 or more vectors by generating the following gnuplot commands:
        // $data << EOD
        // <data values>
        // EOD
        // plot $data using 1:2, $data using 1:3, ...
        template <typename T, typename... Ts>
        void plot(const T& x, Ts... values) {
            db("data", x, values...);
            (*this)("plot ",false);
            doplot<2,Ts...>("$data", values...); // note we peeled off the x values
        }
        
        ~Gnuplot() {
            if (f) {
                fclose(f);
                waitpid(pid,0,0);
            }
            if (ftee) fclose(ftee);
        }
        
        static void test() {
            std::vector<double> x = {1.0,2.0,3.0};
            std::vector<double> y = {-1.0,-2.0,3.0};
            std::vector<double> z = {10.0,11.0,12.0};
            {
                Gnuplot g("set style data lp; set grid");
                g.plot(x,y,z);
            }
            {
                Gnuplot g;
                //g("set term png");
                //g("set output \"test.png\"");
                g.db("xyz",x,y,z);
                g("plot $xyz using 1:2 with linespoints");
            }
            {
                // use x11 to work around bug (https://sourceforge.net/p/gnuplot/bugs/2634/) ... wxt temporarily needs GDK_BACKEND=x11
                Gnuplot g("set term x11; set xrange [-10:10]; set yrange [0:1]; set grid; set style data l", "test3.gnuplot");
                size_t npts = 100;
                std::vector<double> x(npts), y(npts);
                for (size_t i = 0; i<npts; ++i) x[i] = -10.0 + 20.0 * i / (npts-1);
                for (int step=0; step<40; step++) {
                    double phase = step*0.3;
                    for (size_t i = 0; i<npts; ++i) y[i] = 0.5 + 0.5 * std::sin(x[i]+phase);
                    char buf[256];
                    snprintf(buf,sizeof(buf),"set title \"step %d\"",step);
                    g(buf);
                    g.plot(x,y);
                    usleep(100000);
                }
            }
        }
    };
}
#endif
