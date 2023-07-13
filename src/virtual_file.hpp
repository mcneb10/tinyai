#ifndef _VIRTUAL_FILE_HPP_
#define _VIRTUAL_FILE_HPP_

//#include <cstdio>
#include <stdint.h>
#include <string>
#include <cstdarg>

namespace std {

class virtual_file {
		private:
			uint32_t offset = 0;
			string data;
            string data_backup;
		public:
            string name;
			virtual_file(string n, string d) {name=n;data=d;data_backup=d;}
            virtual_file(string n) {name=n;data = "";data_backup="";}
            int scanf(char* fmt, ...) {
                va_list args;
                va_start(args, fmt);

                int r = vsscanf(data.c_str(), fmt, args);
                offset += r;
                va_end(args);
                seek(offset);
                return r;
            }
            void rewind() {seek(0);}
            void seek(uint32_t o) {
                offset = o;
                data = data_backup.substr(offset);
            }
            string get_str() {
                string r;
                while(data.at(0) != '\n') {
                    r += data.at(0);
                    seek(offset+1);
                }
                // Get rid of the new line
                seek(offset+1);
                return r;
            }
            void concat(string s) {
                data_backup += s;
                data = data_backup.substr(offset);
            }
            string dump() {
                return data_backup;
            }
};

}
#endif