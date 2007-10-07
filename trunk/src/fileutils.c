#include <stdio.h>
#include <stdarg.h>

/****************************************************************/
void copyfile(FILE *infp, FILE *outfp)
{
	int c;
	while((c=getc(infp)) != EOF)  {
		putc(c, outfp);
	}
}

/****************************************************************/
// printf to the beginning of the file.
//
// this function takes O(c) disk IOs where c is the size of file fp
//
/****************************************************************/
void finsertf(FILE *fp, char *format, ...)
{
	FILE *tmpf;
	va_list args;

	tmpf = tmpfile(); // creates a temp file

	rewind(fp);
	copyfile(fp, tmpf); // copy from fp to tmpfile

	// insert stuff
	rewind(fp);

	// fprintf
	va_start(args, format);
	vfprintf(fp, format, args);
	va_end(args);

	rewind(tmpf);
	copyfile(tmpf, fp); // copy back
	// no need to clear the old content because EVERYTHING
	// gets overridden.

	close(tmpf); // tmp file gets automatically deleted
}

