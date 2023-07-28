#ifndef _VIRTUAL_FILE_HPP_
#define _VIRTUAL_FILE_HPP_

#include <cstdio>
#include <stdint.h>
#include <string>
#include <cstring>
#include <cstdarg>
#include <cmath>

namespace std
{

    class virtual_file
    {
    private:
        uint32_t offset = 0;
        string data;
        string data_backup;

    public:
        string name;
        virtual_file(string n, string d)
        {
            name = n;
            data = d;
            data_backup = d;
        }
        virtual_file(string n)
        {
            name = n;
            data = "";
            data_backup = "";
        }

        // This is cursed
        // Only use it with ints, uints, or floats otherwise you'll die
        int scanf(const char *fmt, ...)
        {
            va_list args;
            va_list args2;
            va_start(args, fmt);
            va_copy(args2, args);

            int field_count = vsscanf(data.c_str(), fmt, args);
            int i = field_count; //- 1;
            int x;
            // Fuck this
            for (unsigned j = 0; j < strlen(fmt) - 1; j++)
            {
                if (fmt[j] == '%' && fmt[j + 1] != '%')
                {
                    // CURSED
                    switch (fmt[j + 1])
                    {
                    case 'u':
                    {
                        unsigned int arg = *va_arg(args2, unsigned int *);
                        x = arg;
                        i += (x < 10 ? 1 : (x < 100 ? 2 : (x < 1000 ? 3 : (x < 10000 ? 4 : (x < 100000 ? 5 : (x < 1000000 ? 6 : (x < 10000000 ? 7 : (x < 100000000 ? 8 : (x < 1000000000 ? 9 : 10)))))))));
                    }
                    break;
                    case 'i':
                    {
                        int arg = *va_arg(args2, int *);
                        if (arg < 0)
                            i++;
                        x = abs(arg);
                        i += (x < 10 ? 1 : (x < 100 ? 2 : (x < 1000 ? 3 : (x < 10000 ? 4 : (x < 100000 ? 5 : (x < 1000000 ? 6 : (x < 10000000 ? 7 : (x < 100000000 ? 8 : (x < 1000000000 ? 9 : 10)))))))));
                    }
                    break;
                    case 'l':
                        if (fmt[j + 2] == 'f')
                        {
                            double arg = *va_arg(args2, double *);
                            char *double_str = new char[20];
                            sprintf(double_str, "%lf", arg);
                            i += strlen(double_str);
                        }
                        break;
                    default:
                        printf("Unrecognized format in virtual_file.scanf: %%%c\n", fmt[j + 1]);
                        exit(-1);
                    }
                }
            }

            seek(offset + i);
            va_end(args);
            return i;
        }

        void rewind() { seek(0); }
        void seek(uint32_t o)
        {
            offset = o;
            data = data_backup.substr(offset);
        }
        string get_str()
        {
            string r;
            while (data.at(0) != '\n')
            {
                r += data.at(0);
                seek(offset + 1);
            }
            // Get rid of the new line
            seek(offset + 1);
            return r;
        }
        void skip_space_and_newline()
        {
            while (data.at(0) == ' ' || data.at(0) == '\t' || data.at(0) == '\r' || data.at(0) == '\n')
            {
                seek(offset + 1);
            }
        }
        void concat(string s)
        {
            data_backup += s;
            data = data_backup.substr(offset);
        }
        string dump()
        {
            return data_backup;
        }
    };

}
#endif