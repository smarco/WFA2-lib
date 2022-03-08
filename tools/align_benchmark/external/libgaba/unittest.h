
/**
 * @file unittest.h
 *
 * @brief single-header unittesting framework for pure C99 codes
 *
 * @author Hajime Suzuki
 * @date 2016/3/1
 * @license MIT
 *
 * @detail
 * latest version is found at https://github.com/ocxtal/unittest.h.
 * see README.md for the details.
 */
#pragma once	/* instead of include guard */

#ifndef UNITTEST
#define UNITTEST 				1
#endif

#ifndef UNITTEST_ALIAS_MAIN
#define UNITTEST_ALIAS_MAIN		0
#endif

/* for compatibility with -std=c99 (2016/4/26 by Hajime Suzuki) */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

/* end */

#include <alloca.h>
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stddef.h>				/* offsetof */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef UNITTEST_UNIQUE_ID
#define UNITTEST_UNIQUE_ID		0
#endif

/**
 * from kvec.h in klib (https://github.com/attractivechaos/klib)
 * a little bit modified from the original
 */
#define utkv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#define UNITTEST_KV_INIT 			( 64 )

/**
 * output coloring
 */
#define UT_RED				"\x1b[31m"
#define UT_GREEN			"\x1b[32m"
#define UT_YELLOW			"\x1b[33m"
#define UT_BLUE				"\x1b[34m"
#define UT_MAGENTA			"\x1b[35m"
#define UT_CYAN				"\x1b[36m"
#define UT_WHITE			"\x1b[37m"
#define UT_DEFAULT_COLOR	"\x1b[39m"
#define ut_color(c, x)		c x UT_DEFAULT_COLOR

/* static assertion macros */
#define ut_sa_cat_intl(x, y)	x##y
#define ut_sa_cat(x, y)			ut_sa_cat_intl(x, y)
#define ut_static_assert(expr)	typedef char ut_sa_cat(ut_st, __LINE__)[(expr) ? 1 : -1]

/**
 * basic vectors (utkv_*)
 */
#define utkvec_t(type)    struct { uint64_t n, m; type *a; }
#define utkv_init(v)      ( (v).n = 0, (v).m = UNITTEST_KV_INIT, (v).a = (__typeof__((v).a))calloc((v).m, sizeof(*(v).a)) )
#define utkv_destroy(v)   { free((v).a); (v).a = NULL; }
// #define utkv_A(v, i)      ( (v).a[(i)] )
#define utkv_pop(v)       ( (v).a[--(v).n] )
#define utkv_size(v)      ( (v).n )
#define utkv_max(v)       ( (v).m )

#define utkv_clear(v)		( (v).n = 0 )
#define utkv_resize(v, s) ( \
	(v).m = (s), (v).a = (__typeof__((v).a))realloc((v).a, sizeof(*(v).a) * (v).m) )

#define utkv_reserve(v, s) ( \
	(v).m > (s) ? 0 : ((v).m = (s), (v).a = (__typeof__((v).a))realloc((v).a, sizeof(*(v).a) * (v).m), 0) )

#define utkv_copy(v1, v0) do {								\
		if ((v1).m < (v0).n) utkv_resize(v1, (v0).n);			\
		(v1).n = (v0).n;									\
		memcpy((v1).a, (v0).a, sizeof(*(v).a) * (v0).n);	\
	} while (0)												\

#define utkv_push(v, x) do {									\
		if ((v).n == (v).m) {								\
			(v).m = (v).m * 2;								\
			(v).a = (__typeof__((v).a))realloc((v).a, sizeof(*(v).a) * (v).m);	\
		}													\
		(v).a[(v).n++] = (x);								\
	} while (0)

#define utkv_pushp(v) ( \
	((v).n == (v).m) ?									\
	((v).m = (v).m * 2,									\
	 (v).a = realloc((v).a, sizeof(*(v).a) * (v).m), 0)	\
	: 0), ( (v).a + ((v).n++) )

#define utkv_a(v, i) ( \
	((v).m <= (size_t)(i) ? \
	((v).m = (v).n = (i) + 1, utkv_roundup32((v).m), \
	 (v).a = realloc((v).a, sizeof(*(v).a) * (v).m), 0) \
	: (v).n <= (size_t)(i) ? (v).n = (i) + 1 \
	: 0), (v).a[(i)])

/** bound-unchecked accessor */
#define utkv_at(v, i) ( (v).a[(i)] )
#define utkv_ptr(v)  ( (v).a )

/**
 * forward declaration of the main function
 */
int main(int argc, char *argv[]);

/**
 * @struct ut_result_s
 */
struct ut_result_s {
	int64_t cnt;
	int64_t succ;
	int64_t fail;
};

struct ut_global_config_s;
struct ut_group_config_s;
struct ut_s;
struct ut_result_s;
struct ut_printer_s {
	/* printers */
	void (*failed)(
		struct ut_s const *info,
		struct ut_global_config_s const *gconf,
		struct ut_group_config_s const *config,
		int64_t line,
		char const *func,
		char const *expr,
		char const *fmt,
		...);

	void (*result)(
		struct ut_global_config_s const *gconf,
		struct ut_group_config_s const *config,
		struct ut_result_s const *result,
		int64_t file_cnt);
};

/**
 * @struct ut_global_config_s
 */
struct ut_global_config_s {
	FILE *fp;
	struct ut_printer_s printer;
	uint64_t threads;
};

/**
 * @struct ut_group_config_s
 *
 * @brief unittest scope config
 */
struct ut_group_config_s {
	/* internal use */
	char const *file;
	int64_t unique_id;
	uint64_t line;
	int64_t exec;

	uint64_t _pad[4];

	/* dependency resolution */
	char const *name;
	char const *depends_on[16];

	/* environment setup and cleanup */
	void *(*init)(void *params);
	void (*clean)(void *context);
	void *params;
};

/**
 * @struct ut_s
 *
 * @brief unittest function config
 */
struct ut_s {
	/* for internal use */
	char const *file;
	int64_t unique_id;
	uint64_t line;
	int64_t exec;

	/* per-function config */
	void (*fn)(
		void *ctx,
		void *gctx,
		struct ut_s *info,
		struct ut_global_config_s const *ut_gconf,
		struct ut_group_config_s const *config);

	/* internal use (cont'd) */
	uint64_t index;
	uint64_t succ, fail;		/* result: 1 when succeeded, 2 when failed */

	/* dependency resolution */
	char const *name;
	char const *depends_on[16];

	/* environment setup and cleanup */
	void *(*init)(void *params);
	void (*clean)(void *context);
	void *params;
};

/* the two structs must be castable */
ut_static_assert(offsetof(struct ut_group_config_s, file) == offsetof(struct ut_s, file));
ut_static_assert(offsetof(struct ut_group_config_s, unique_id) == offsetof(struct ut_s, unique_id));
ut_static_assert(offsetof(struct ut_group_config_s, line) == offsetof(struct ut_s, line));
ut_static_assert(offsetof(struct ut_group_config_s, exec) == offsetof(struct ut_s, exec));
ut_static_assert(offsetof(struct ut_group_config_s, name) == offsetof(struct ut_s, name));
ut_static_assert(offsetof(struct ut_group_config_s, depends_on) == offsetof(struct ut_s, depends_on));
ut_static_assert(offsetof(struct ut_group_config_s, init) == offsetof(struct ut_s, init));
ut_static_assert(offsetof(struct ut_group_config_s, clean) == offsetof(struct ut_s, clean));
ut_static_assert(offsetof(struct ut_group_config_s, params) == offsetof(struct ut_s, params));


/**
 * @macro ut_unused
 * @brief declare the variable is unused in the function
 */
#define ut_unused(x)						(void)(x);

/**
 * @macro ut_null_replace
 * @brief replace pointer with "(null)"
 */
#define ut_null_replace(ptr, str)			( (ptr) == NULL ? (str) : (ptr) )

/**
 * @macro ut_build_name
 *
 * @brief an utility macro to make unique name
 */
#define ut_join_name(a, b, c)				a##b##_##c
#define ut_build_name(prefix, num, id)	ut_join_name(prefix, num, id)

/**
 * @macro unittest
 *
 * @brief instanciate a unittest object
 */
#define UNITTEST_ARG_DECL \
	void *ctx, \
	void *gctx, \
	struct ut_s *ut_info, \
	struct ut_global_config_s const *ut_gconf, \
	struct ut_group_config_s const *ut_config
#define UNITTEST_ARG_LIST 	ctx, gctx, ut_info, ut_gconf, ut_config
#if UNITTEST != 0
#define unittest(...) \
	static void ut_build_name(ut_body_, UNITTEST_UNIQUE_ID, __LINE__)(UNITTEST_ARG_DECL); \
	static struct ut_s const ut_build_name(ut_info_, UNITTEST_UNIQUE_ID, __LINE__) = { \
		/* file, unique_id, line, exec, fn, name, depends_on */ \
		__FILE__, UNITTEST_UNIQUE_ID, __LINE__, 1, ut_build_name(ut_body_, UNITTEST_UNIQUE_ID, __LINE__), \
		0, 0, 0, __VA_ARGS__ \
	}; \
	struct ut_s ut_build_name(ut_get_info_, UNITTEST_UNIQUE_ID, __LINE__)(void) \
	{ \
		return(ut_build_name(ut_info_, UNITTEST_UNIQUE_ID, __LINE__)); \
	} \
	static void ut_build_name(ut_body_, UNITTEST_UNIQUE_ID, __LINE__)(UNITTEST_ARG_DECL)

#else	/* UNITTEST != 0 */

#define unittest(...) \
	static void ut_build_name(ut_body_, UNITTEST_UNIQUE_ID, __LINE__)(UNITTEST_ARG_DECL)


#endif	/* UNITTEST != 0 */

/**
 * @macro unittest_config
 *
 * @brief scope configuration
 */
#if UNITTEST != 0
#define unittest_config(...) \
	static struct ut_group_config_s const ut_build_name(ut_config_, UNITTEST_UNIQUE_ID, __LINE__) = { \
		.file = __FILE__, \
		.line = __LINE__, \
		.unique_id = UNITTEST_UNIQUE_ID, \
		__VA_ARGS__ \
	}; \
	struct ut_group_config_s ut_build_name(ut_get_config_, UNITTEST_UNIQUE_ID, 0)(void) \
	{ \
		return(ut_build_name(ut_config_, UNITTEST_UNIQUE_ID, __LINE__)); \
	}

#else	/* UNITTEST != 0 */

#define unittest_config(...)

#endif	/* UNITTEST != 0 */

/* assertion failed message printers */
static
void ut_print_assertion_failed(
	struct ut_s const *info,
	struct ut_global_config_s const *gconf,
	struct ut_group_config_s const *config,
	int64_t line,
	char const *func,
	char const *expr,
	char const *fmt,
	...)
{
	ut_unused(func);

	va_list l;
	va_start(l, fmt);

	fprintf(gconf->fp,
		ut_color(UT_YELLOW, "assertion failed") ": [%s] %s:" ut_color(UT_BLUE, "%" PRId64 "") " (%s) `" ut_color(UT_MAGENTA, "%s") "'",
		ut_null_replace(config->name, "no name"),
		ut_null_replace(info->file, "(unknown filename)"),
		line,
		ut_null_replace(info->name, "no name"),
		expr);
	if(strlen(fmt) != 0) {
		fprintf(gconf->fp, ", ");
		vfprintf(gconf->fp, fmt, l);
	}
	fprintf(gconf->fp, "\n");
	va_end(l);
	return;
}

static
void ut_print_assertion_failed_json(
	struct ut_s const *info,
	struct ut_global_config_s const *gconf,
	struct ut_group_config_s const *config,
	int64_t line,
	char const *func,
	char const *expr,
	char const *fmt,
	...)
{
	ut_unused(config);
	ut_unused(func);

	va_list l;
	va_start(l, fmt);

	fprintf(gconf->fp, "{ \"tag\": \"fail\", ");
	if(config->name != NULL) {
		fprintf(gconf->fp, "\"group\": \"%s\", ", config->name);
	}
	if(config->file != NULL) {
		fprintf(gconf->fp, "\"filename\": \"%s\", ", config->file);
	}

	fprintf(gconf->fp, "\"line\": %" PRId64 ", ", line);
	if(info->name != NULL) {
		fprintf(gconf->fp, "\"name\": \"%s\", ", info->name);
	}
	fprintf(gconf->fp, "\"expr\": \"%s\", ", expr);
	if(strlen(fmt) != 0) {
		fprintf(gconf->fp, "\"debugprint\": \"");
		vfprintf(gconf->fp, fmt, l);			/* FIXME: escape */
		fprintf(gconf->fp, "\", ");
	}
	fprintf(gconf->fp, "},\n");
	va_end(l);
	return;
}

/* summary printers */
static
void ut_print_results(
	struct ut_global_config_s const *gconf,
	struct ut_group_config_s const *config,
	struct ut_result_s const *result,
	int64_t file_cnt)
{
	int64_t cnt = 0;
	int64_t succ = 0;
	int64_t fail = 0;

	for(int64_t i = 0, j = 0; i < file_cnt; i++) {
		if(config[i].exec == 0) { continue; }

		fprintf(gconf->fp, "%sGroup %s: %" PRId64 " succeeded, %" PRId64 " failed in total %" PRId64 " assertions in %" PRId64 " tests.%s\n",
			(result[j].fail == 0) ? UT_GREEN : UT_RED,
			ut_null_replace(config[i].name, "(no name)"),
			result[j].succ,
			result[j].fail,
			result[j].succ + result[j].fail,
			result[j].cnt,
			UT_DEFAULT_COLOR);
		
		cnt += result[j].cnt;
		succ += result[j].succ;
		fail += result[j].fail;
		j++;
	}

	fprintf(gconf->fp, "%sSummary: %" PRId64 " succeeded, %" PRId64 " failed in total %" PRId64 " assertions in %" PRId64 " tests.%s\n",
		(fail == 0) ? UT_GREEN : UT_RED,
		succ, fail, succ + fail, cnt,
		UT_DEFAULT_COLOR);
	return;
}

static
void ut_print_results_json(
	struct ut_global_config_s const *gconf,
	struct ut_group_config_s const *config,
	struct ut_result_s const *result,
	int64_t file_cnt)
{
	int64_t cnt = 0;
	int64_t succ = 0;
	int64_t fail = 0;

	fprintf(gconf->fp, "{ \"tag\": \"results\", [ ");

	for(int64_t i = 0, j = 0; i < file_cnt; i++) {
		if(config[i].exec == 0) { continue; }
		if(config->name != NULL) {
			fprintf(gconf->fp, "{ \"group\": \"%s\", ", config->name);
		}
		if(config->file != NULL) {
			fprintf(gconf->fp, "\"filename\": \"%s\", ", config->file);
		}
		fprintf(gconf->fp, "\"succeeded\": %" PRId64 ", ", result[j].succ);
		fprintf(gconf->fp, "\"failed\": %" PRId64 ", ", result[j].fail);
		fprintf(gconf->fp, "\"assertioncount\": %" PRId64 ", ", result[j].succ + result[j].fail);
		fprintf(gconf->fp, "\"testcount\": %" PRId64 ", ", result[j].cnt);
		fprintf(gconf->fp, "}, ");
		
		cnt += result[j].cnt;
		succ += result[j].succ;
		fail += result[j].fail;
		j++;
	}
	fprintf(gconf->fp, "] },\n");


	fprintf(gconf->fp, "{ \"tag\": \"summary\", ");
	fprintf(gconf->fp, "\"succeeded\": %" PRId64 ", ", succ);
	fprintf(gconf->fp, "\"failed\": %" PRId64 ", ", fail);
	fprintf(gconf->fp, "\"assertioncount\": %" PRId64 ", ", succ + fail);
	fprintf(gconf->fp, "\"testcount\": %" PRId64 ", ", cnt);
	fprintf(gconf->fp, "},\n");

	return;
}

static
struct ut_printer_s ut_default_printer = {
	.failed = ut_print_assertion_failed,
	.result = ut_print_results
};

static
struct ut_printer_s ut_json_printer = {
	.failed = ut_print_assertion_failed_json,
	.result = ut_print_results_json
};

/**
 * memory dump macro
 */
#define ut_dump(ptr, len) ({ \
	uint64_t _size = (((len) + 15) / 16 + 1) * \
		(strlen("0x0123456789abcdef:") + 16 * strlen(" 00a") + strlen("  \n+ margin")) \
		+ strlen(#ptr) + strlen("\n`' len: 100000000"); \
	uint8_t *_ptr = (uint8_t *)(ptr); \
	char *_str = alloca(_size); \
	char *_s = _str; \
	/* make header */ \
	_s += sprintf(_s, "\n`%s' len: %" PRIu64 "\n", #ptr, (uint64_t)len); \
	_s += sprintf(_s, "                   "); \
	for(uint64_t i = 0; i < 16; i++) { \
		_s += sprintf(_s, " %02x", (uint8_t)i); \
	} \
	_s += sprintf(_s, "\n"); \
	for(uint64_t i = 0; i < ((len) + 15) / 16; i++) { \
		_s += sprintf(_s, "0x%016" PRIx64 ":", (uint64_t)_ptr); \
		for(uint64_t j = 0; j < 16; j++) { \
			_s += sprintf(_s, " %02x", (uint8_t)_ptr[j]); \
		} \
		_s += sprintf(_s, "  "); \
		for(uint64_t j = 0; j < 16; j++) { \
			_s += sprintf(_s, "%c", isprint(_ptr[j]) ? _ptr[j] : ' '); \
		} \
		_s += sprintf(_s, "\n"); \
		_ptr += 16; \
	} \
	(char const *)_str; \
})
#ifndef dump
// #define dump 			ut_dump
#endif

/**
 * assertion macro
 */
#define ut_assert(expr, ...) { \
	if(expr) { \
		ut_info->succ++; \
	} else { \
		ut_info->fail++; \
		/* dump debug information */ \
		ut_gconf->printer.failed(ut_info, ut_gconf, ut_config, __LINE__, __func__, #expr, "" __VA_ARGS__); \
	} \
}
#ifndef assert
#define assert			ut_assert
#endif

/**
 * @struct ut_nm_result_s
 * @brief parsed result container
 */
struct ut_nm_result_s {
	void *ptr;
	char type;
	char name[255];
};

static inline
int ut_strcmp(
	char const *a,
	char const *b)
{
	/* if both are NULL */
	if(a == NULL && b == NULL) {
		return(0);
	}

	if(b == NULL) {
		return(1);
	}
	if(a == NULL) {
		return(-1);
	}
	return(strcmp(a, b));
}

static inline
char *ut_build_nm_cmd(
	char const *filename)
{
	uint64_t const filename_len_limit = 1024;
	char const *cmd_base = "nm ";

	/* check the length of the filename */
	if(strlen(filename) > filename_len_limit) {
		return(NULL);
	}

	/* cat name */
	char *cmd = (char *)malloc(strlen(cmd_base) + strlen(filename) + 1);
	strcpy(cmd, cmd_base);
	strcat(cmd, filename);

	return(cmd);
}

static inline
char *ut_dump_file(
	FILE *fp)
{
	int c;
	utkvec_t(char) buf;

	utkv_init(buf);
	while((c = getc(fp)) != EOF) {
		utkv_push(buf, c);
	}

	/* push terminator */
	utkv_push(buf, '\0');
	return(utkv_ptr(buf));
}

static inline
char *ut_dump_nm_output(
	char const *filename)
{
	char *cmd = NULL;
	FILE *fp = NULL;
	char *res = NULL;

	/* build command */
	if((cmd = ut_build_nm_cmd(filename)) == NULL) {
		goto _ut_nm_error_handler;
	}

	/* open */
	if((fp = popen(cmd, "r")) == NULL) {
		fprintf(stderr, ut_color(UT_RED, "ERROR") ": failed to open pipe.\n");
		goto _ut_nm_error_handler;
	}

	/* dump */
	if((res = ut_dump_file(fp)) == NULL) {
		fprintf(stderr, ut_color(UT_RED, "ERROR") ": failed to read nm output.\n");
		goto _ut_nm_error_handler;
	}

	/* close file */
	if(pclose(fp) != 0) {
		fprintf(stderr, ut_color(UT_RED, "ERROR") ": failed to close pipe.\n");
		goto _ut_nm_error_handler;
	}
	free(cmd); cmd = NULL;
	return(res);

_ut_nm_error_handler:
	if(cmd != NULL) { free(cmd); }
	if(fp != NULL) { pclose(fp); }
	if(res != NULL) { free(res); }
	return(NULL);
}

static inline
struct ut_nm_result_s *ut_parse_nm_output(
	char const *str)
{
	char const *p = str;
	utkvec_t(struct ut_nm_result_s) buf;

	utkv_init(buf);
	while(*p != '\0') {
		struct ut_nm_result_s r;

		/* check the sanity of the line */
		char const *sp = p;
		while(*sp != '\r' && *sp != '\n' && *sp != '\0') { sp++; }
		// printf("%p, %p, %" PRId64 "\n", p, sp, (int64_t)(sp - p));
		if((int64_t)(sp - p) < 1) { break; }

		/* if the first character of the line is space, pass the PTR state */
		if(!isspace(*p)) {
			char *np;
			r.ptr = (void *)((uint64_t)strtoll(p, &np, 16));

			/* advance pointer */
			p = np;

			// printf("ptr field found: %p\n", r.ptr);
		} else {
			r.ptr = NULL;

			// printf("ptr field not found\n");
		}

		/* parse type */
		while(isspace(*p)) { p++; }
		r.type = *p++;

		// printf("type found %c\n", r.type);

		/* parse name */
		while(isspace(*p)) { p++; }
		if(*p == '_') { p++; }
		int name_idx = 0;
		while(*p != '\n' && name_idx < 254) {
			r.name[name_idx++] = *p++;
		}
		r.name[name_idx] = '\0';

		// printf("name found %s\n", r.name);

		/* adjust the pointer to the head of the next line */
		while(*p == '\r' || *p == '\n') { p++; }

		/* push */
		utkv_push(buf, r);
	}

	/* push terminator */
	utkv_push(buf, (struct ut_nm_result_s){ 0 });
	return(utkv_ptr(buf));
}

static inline
struct ut_nm_result_s *ut_nm(
	char const *filename)
{
	char const *str = NULL;
	struct ut_nm_result_s *res = NULL;

	if((str = ut_dump_nm_output(filename)) == NULL) {
		return(NULL);
	}
	res = ut_parse_nm_output(str);
	
	free((void *)str);
	return(res);
}

static inline
int ut_startswith(char const *str, char const *prefix)
{
	while(*str != '\0' && *prefix != '\0') {
		if(*str++ != *prefix++) { return(1); }
	}
	return(*prefix == '\0' ? 0 : 1);
}

static inline
struct ut_s *ut_get_unittest(
	struct ut_nm_result_s const *res)
{
	/* extract offset */
	uintptr_t offset = (uintptr_t)-1LL;
	struct ut_nm_result_s const *r = res;
	while(r->type != (char)0) {
		if(ut_strcmp("main", r->name) == 0) {
			offset = (uintptr_t)main - (uintptr_t)r->ptr;
		}
		r++;
	}

	if(offset == (uintptr_t)-1LL) {
		return(NULL);
	}

	#define ut_get_info_call_func(_ptr, _offset) ( \
		(struct ut_s (*)(void))((uintptr_t)(_ptr) + (uintptr_t)(_offset)) \
	)

	/* get info */
	utkvec_t(struct ut_s) buf;
	
	r = res;
	utkv_init(buf);
	while(r->type != (char)0) {
		if(ut_startswith(r->name, "ut_get_info_") == 0) {
			struct ut_s i = ut_get_info_call_func(r->ptr, offset)();
			utkv_push(buf, i);
		}
		r++;
	}

	/* push terminator */
	utkv_push(buf, (struct ut_s){ 0 });
	return(utkv_ptr(buf));
}

static inline
struct ut_group_config_s *ut_get_ut_config(
	struct ut_nm_result_s const *res)
{
	/* extract offset */
	uintptr_t offset = (uintptr_t)-1LL;
	struct ut_nm_result_s const *r = res;
	while(r->type != (char)0) {
		if(ut_strcmp("main", r->name) == 0) {
			offset = (uintptr_t)main - (uintptr_t)r->ptr;
		}
		r++;
	}

	if(offset == (uintptr_t)-1LL) {
		return(NULL);
	}

	#define ut_get_config_call_func(_ptr, _offset) ( \
		(struct ut_group_config_s (*)(void))((uintptr_t)(_ptr) + (uintptr_t)(_offset)) \
	)

	/* get info */
	utkvec_t(struct ut_group_config_s) buf;
	
	r = res;
	utkv_init(buf);
	while(r->type != (char)0) {
		if(ut_startswith(r->name, "ut_get_config_") == 0) {
			struct ut_group_config_s i = ut_get_config_call_func(r->ptr, offset)();
			utkv_push(buf, i);
		}
		r++;
	}

	/* push terminator */
	utkv_push(buf, (struct ut_group_config_s){ 0 });
	return(utkv_ptr(buf));
}

static inline
void ut_dump_test(
	struct ut_s const *test)
{
	struct ut_s const *t = test;
	while(t->file != NULL) {
		printf("%s, %" PRIu64 ", %" PRId64 ", %s, %s, %p, %p\n",
			t->file,
			t->line,
			t->unique_id,
			t->name,
			t->depends_on[0],
			(void *)t->init,
			(void *)t->clean);
		t++;
	}
	return;
}

static inline
void ut_dump_config(
	struct ut_group_config_s const *config)
{
	struct ut_group_config_s const *c = config;
	while(c->file != NULL) {
		printf("%s, %" PRId64 ", %s, %s, %p, %p\n",
			c->file,
			c->unique_id,
			c->name,
			c->depends_on[0],
			(void *)c->init,
			(void *)c->clean);
		c++;
	}
	return;
}

static
int ut_compare(
	void const *_a,
	void const *_b)
{
	struct ut_s const *a = (struct ut_s const *)_a;
	struct ut_s const *b = (struct ut_s const *)_b;

	int64_t comp_res = 0;
	/* first sort by file name */
	if((comp_res = ut_strcmp(a->file, b->file)) != 0) {
		return(comp_res);
	}

	#define sat(a)	( ((a) > INT32_MAX) ? INT32_MAX : (((a) < INT32_MIN) ? INT32_MIN : (a)) )

	if((comp_res = (a->unique_id - b->unique_id)) != 0) {
		return(sat(comp_res));
	}

	#undef sat

	/* third sort by name */
	if((comp_res = ut_strcmp(a->name, b->name)) != 0) {
		return(comp_res);
	}

	/* last, sort by line number */
	return((int)(a->line - b->line));
}

static inline
int ut_match(
	void const *_a,
	void const *_b)
{
	struct ut_s const *a = (struct ut_s *)_a;
	struct ut_s const *b = (struct ut_s *)_b;
	return((ut_strcmp(a->file, b->file) == 0 && a->unique_id == b->unique_id) ? 0 : 1);
}

static inline
uint64_t ut_get_total_test_count(
	struct ut_s const *test)
{
	uint64_t cnt = 0;
	struct ut_s const *t = test;
	while(t->file != NULL) { t++; cnt++; }
	return(cnt);
}

static inline
uint64_t ut_get_total_config_count(
	struct ut_group_config_s const *config)
{
	uint64_t cnt = 0;
	struct ut_group_config_s const *c = config;
	while(c->file != NULL) { c++; cnt++; }
	return(cnt);
}

static inline
uint64_t ut_get_total_file_count(
	struct ut_s const *sorted_test)
{
	uint64_t cnt = 1;
	struct ut_s const *t = sorted_test;

	if(t++->file == NULL) {
		return(0);
	}

	while(t->file != NULL) {
		// printf("comp %s and %s\n", (t - 1)->file, t->file);
		if(ut_match((void *)(t - 1), (void *)t) != 0) {
			cnt++;
		}
		t++;
	}
	return(cnt);
}

static inline
void ut_sort(
	struct ut_s *test,
	struct ut_group_config_s *config)
{
	// ut_dump_test(test);
	qsort(test,
		ut_get_total_test_count(test),
		sizeof(struct ut_s),
		ut_compare);
	qsort(config,
		ut_get_total_config_count(config),
		sizeof(struct ut_group_config_s),
		ut_compare);
	// ut_dump_test(test);
	// printf("%" PRId64 "\n", ut_get_total_file_count(test));
	return;
}

static inline
uint64_t *ut_build_file_index(
	struct ut_s const *sorted_test)
{
	uint64_t cnt = 1;
	utkvec_t(uint64_t) idx;
	struct ut_s const *t = sorted_test;

	utkv_init(idx);
	if(t++->file == NULL) {
		utkv_push(idx, 0);
		utkv_push(idx, -1);
		return(utkv_ptr(idx));
	}

	utkv_push(idx, 0);
	while(t->file != NULL) {
		// printf("comp %s and %s\n", (t - 1)->name, t->name);
		if(ut_match((void *)(t - 1), (void *)t) != 0) {
			utkv_push(idx, cnt);
		}
		cnt++;
		t++;
	}

	utkv_push(idx, cnt);
	utkv_push(idx, -1);
	return(utkv_ptr(idx));
}

static inline
struct ut_group_config_s *ut_compensate_config(
	struct ut_s *sorted_test,
	struct ut_group_config_s *sorted_config)
{
	// printf("compensate config\n");
	// ut_dump_test(sorted_test);
	// ut_dump_config(sorted_config);

	uint64_t *file_idx = ut_build_file_index(sorted_test);
	utkvec_t(struct ut_group_config_s) compd_config;
	utkv_init(compd_config);

	int64_t i = 0;
	// printf("%" PRId64 "\n", file_idx[i]);
	struct ut_s const *t = &sorted_test[file_idx[i]];
	struct ut_group_config_s const *c = sorted_config;

	while(c->file != NULL && t->file != NULL) {
		while(ut_match((void *)t, (void *)c) != 0) {
			// printf("filename does not match %s, %s\n", c->file, t->file);
			utkv_push(compd_config, (struct ut_group_config_s){ .file = t->file });
			t = &sorted_test[file_idx[++i]];
		}
		// printf("matched %s, %s\n", c->file, t->file);
		utkv_push(compd_config, *c);
		c++;
		t = &sorted_test[file_idx[++i]];
	}

	free(file_idx);
	return(utkv_ptr(compd_config));
}

static inline
int ut_toposort_by_tag(
	struct ut_s *sorted_test,
	uint64_t test_cnt)
{
	/* init dag */
	utkvec_t(utkvec_t(uint64_t)) dag;
	utkv_init(dag);
	utkv_reserve(dag, test_cnt);
	for(uint64_t i = 0; i < test_cnt; i++) {
		utkv_init(utkv_at(dag, i));
	}

	/* build dag */
	for(uint64_t i = 0; i < test_cnt; i++) {
		/* enumerate depends_on */
		char const *const *d = sorted_test[i].depends_on;
		while(*d != NULL) {
			/* enumerate tests */
			for(uint64_t j = 0; j < test_cnt; j++) {
				if(i == j) { continue; }
				if(ut_strcmp(sorted_test[j].name, *d) == 0) {
					// printf("edge found from node %" PRId64 " to node %" PRId64 "\n", j, i);
					utkv_push(utkv_at(dag, i), j);
				}
			}
			d++;
		}
	}

	/* build mark */
	utkvec_t(int8_t) mark;
	utkv_init(mark);
	for(uint64_t i = 0; i < test_cnt; i++) {
		// printf("incoming edge count %" PRId64 " : %" PRId64 "\n", i, utkv_size(utkv_at(dag, i)));
		utkv_push(mark, utkv_size(utkv_at(dag, i)));
	}

	/* sort */
	utkvec_t(struct ut_s) res;
	utkv_init(res);
	for(uint64_t i = 0; i < test_cnt; i++) {
		uint64_t node_id = (uint64_t)-1LL;
		for(uint64_t j = 0; j < test_cnt; j++) {
			/* check if the node has incoming edges or is already pushed */
			if(utkv_at(mark, j) == 0) {
				node_id = j; break;
			}
		}

		if(node_id == (uint64_t)-1LL) {
			fprintf(stderr,
				ut_color(UT_RED, "ERROR") ": detected circular dependency in the tests in `" ut_color(UT_MAGENTA, "%s") "'.\n",
				sorted_test[0].file);
			return(-1);
		}

		/* the node has no incoming edge */
		// printf("the node %" PRId64 " has no incoming edge\n", node_id);
		utkv_push(res, sorted_test[node_id]);

		/* mark pushed */
		utkv_at(mark, node_id) = -1;

		/* delete edges */
		for(uint64_t j = 0; j < test_cnt; j++) {
			if(node_id == j) { continue; }

			for(uint64_t k = 0; k < utkv_size(utkv_at(dag, j)); k++) {
				if(utkv_at(utkv_at(dag, j), k) == node_id) {
					// printf("the node %" PRId64 " had edge from %" PRId64 "\n", j, node_id);
					utkv_at(utkv_at(dag, j), k) = -1;
					utkv_at(mark, j)--;
					// printf("and then node %" PRId64 " has %d incoming edges\n", j, utkv_at(mark, j));
				}
			}
		}
	}

	if(utkv_size(res) != test_cnt) {
		fprintf(stderr,
			ut_color(UT_RED, "ERROR") ": detected circular dependency in the tests in `" ut_color(UT_MAGENTA, "%s") "'.\n",
			sorted_test[0].file);
		return(-1);
	}

	/* write back */
	for(uint64_t i = 0; i < test_cnt; i++) {
		sorted_test[i] = utkv_at(res, i);
	}

	/* cleanup */
	for(uint64_t i = 0; i < test_cnt; i++) {
		utkv_destroy(utkv_at(dag, i));
	}
	utkv_destroy(dag);
	utkv_destroy(mark);
	utkv_destroy(res);

	return(0);
}

static inline
int ut_toposort_by_group(
	struct ut_s *sorted_test,
	uint64_t test_cnt,
	struct ut_group_config_s *sorted_config,
	uint64_t *file_idx,
	uint64_t file_cnt)
{
	/* init dag */
	utkvec_t(utkvec_t(uint64_t)) dag;
	utkv_init(dag);
	utkv_reserve(dag, file_cnt);
	for(uint64_t i = 0; i < file_cnt; i++) {
		utkv_init(utkv_at(dag, i));
	}

	/* build dag */
	for(uint64_t i = 0; i < file_cnt; i++) {
		/* enumerate depends_on */
		char const *const *d = sorted_config[i].depends_on;
		while(*d != NULL) {
			/* enumerate tests */
			for(uint64_t j = 0; j < file_cnt; j++) {
				if(i == j) { continue; }
				if(ut_strcmp(sorted_config[j].name, *d) == 0) {
					// printf("edge found from group %" PRId64 " to group %" PRId64 "\n", j, i);
					utkv_push(utkv_at(dag, i), j);
				}
			}
			d++;
		}
	}

	/* build mark */
	utkvec_t(int8_t) mark;
	utkv_init(mark);
	for(uint64_t i = 0; i < file_cnt; i++) {
		// printf("incoming edge count %" PRId64 " : %" PRId64 "\n", i, utkv_size(utkv_at(dag, i)));
		utkv_push(mark, utkv_size(utkv_at(dag, i)));
	}

	/* sort */
	utkvec_t(struct ut_s) test_buf;
	utkvec_t(struct ut_group_config_s) config_buf;
	utkv_init(test_buf);
	utkv_init(config_buf);
	// utkvec_t(uint64_t) res;
	// utkv_init(res);
	for(uint64_t i = 0; i < file_cnt; i++) {
		uint64_t file_id = (uint64_t)-1LL;
		for(uint64_t j = 0; j < file_cnt; j++) {
			/* check if the node has incoming edges or is already pushed */
			if(utkv_at(mark, j) == 0) {
				file_id = j; break;
			}
		}

		if(file_id == (uint64_t)-1LL) {
			fprintf(stderr, ut_color(UT_RED, "ERROR") ": detected circular dependency between groups.\n");
			return(-1);
		}

		/* the node has no incoming edge */
		// printf("the node %" PRId64 " has no incoming edge\n", file_id);
		// utkv_push(res, file_id);

		utkv_push(config_buf, sorted_config[file_id]);
		for(uint64_t j = file_idx[file_id]; j < file_idx[file_id + 1]; j++) {
			utkv_push(test_buf, sorted_test[j]);
		}

		/* mark pushed */
		utkv_at(mark, file_id) = -1;

		/* delete edges */
		for(uint64_t j = 0; j < file_cnt; j++) {
			if(file_id == j) { continue; }

			for(uint64_t k = 0; k < utkv_size(utkv_at(dag, j)); k++) {
				if(utkv_at(utkv_at(dag, j), k) == file_id) {
					// printf("the node %" PRId64 " had edge from %" PRId64 "\n", j, file_id);
					utkv_at(utkv_at(dag, j), k) = -1;
					utkv_at(mark, j)--;
					// printf("and then node %" PRId64 " has %d incoming edges\n", j, utkv_at(mark, j));
				}
			}
		}
	}

	if(utkv_size(config_buf) != file_cnt) {
		fprintf(stderr, ut_color(UT_RED, "ERROR") ": detected circular dependency between groups.\n");
		return(-1);
	}

	/* write back */
	for(uint64_t i = 0; i < test_cnt; i++) {
		sorted_test[i] = utkv_at(test_buf, i);
	}
	for(uint64_t i = 0; i < file_cnt; i++) {
		sorted_config[i] = utkv_at(config_buf, i);
	}

	/* cleanup */
	for(uint64_t i = 0; i < file_cnt; i++) {
		utkv_destroy(utkv_at(dag, i));
	}
	utkv_destroy(dag);
	utkv_destroy(mark);
	utkv_destroy(test_buf);
	utkv_destroy(config_buf);

	return(0);
}

/**
 * @fn ut_build_short_option_string
 * @brief build short getopt string (i.e. "a:b:vVQ") from an array of struct option
 */
static inline
char *ut_build_short_option_string(struct option const *opts)
{
	char *str = NULL, *ps;
	uint64_t len = 0;
	struct option const *po = NULL;

	for(po = opts; po->name != NULL; po++) { len++; }
	str = ps = (char *)malloc(2 * len + 1);
	for(po = opts; po->name != NULL; po++) {
		*ps++ = (char)po->val;
		if(po->has_arg != no_argument) {
			*ps++ = ':';
		}
	}
	*ps = '\0';
	return(str);
}

/**
 * @fn ut_modify_test_config_mark
 */
static inline
int ut_modify_test_config_mark(
	char const *arg,
	void *_test,
	int64_t cnt)
{
	struct ut_s *test = (struct ut_s *)_test;
	char const *p = arg, *b = arg;
	while(*p != '\0') {
		/* parse with comma */
		while(*p != '\0' && *p != ',') { p++; }

		/* copy string */
		char buf[p - b + 1];
		memcpy(buf, b, p - b);
		buf[p - b] = '\0';

		/* linear search among tests */
		int marked = 0;
		for(int64_t i = 0; i < cnt; i++) {
			if(ut_strcmp(test[i].name, buf) == 0) {
				test[i].exec = 2; marked = 1;
			} else {
				test[i].exec = 0;
			}
		}
		if(marked == 0) {
			fprintf(stderr, ut_color(UT_YELLOW, "Warning") ": group `%s' not found.\n", buf);
		}

		if(*p == '\0') { break; }
		b = ++p;
	}
	return(0);
}

/**
 * @fn ut_modify_test_config_all
 */
static inline
int ut_modify_test_config_all(
	void *_test,
	int64_t cnt)
{
	struct ut_s *test = (struct ut_s *)_test;
	for(int64_t i = 0; i < cnt; i++) {
		test[i].exec = 1;
	}
	return(0);
}

/**
 * @fn ut_print_help
 */
static inline
void ut_print_help(void)
{
	static char const *help_message =
		"\n"
		"  unittest.h -- simple unittesting framework for C99\n"
		"\n"
		"  Options:\n"
		"    -g, --group   [STR,...]  specify group names\n"
		"    -t, --test    [STR,...]  specify test names\n"
		"    -s, --seed    [INT       invoke libc::srand with the seed\n"
		"    -o, --stdout             redirect to stdout\n"
		"    -j, --json               print result in json\n"
		"    -n, --threads [INT]      number of threads\n"
		"    -h, --help               show this message\n"
		"\n"
		"  this is an auto-generated message from unittest.h\n"
		"  the latest source is available at:\n"
		"    https://github.com/ocxtal/unittest.h\n"
		"\n";

	fprintf(stderr, "%s", help_message);
	return;
}

/**
 * @fn ut_modify_test_config
 */
static inline
int ut_modify_test_config(
	int argc,
	char *const argv[],
	struct ut_global_config_s *params,
	struct ut_s *sorted_test,
	int64_t test_cnt,
	struct ut_group_config_s *sorted_config,
	int64_t file_cnt)
{
	static struct option const opts_long[] = {
		{ "group", required_argument, NULL, 'g' },
		{ "test", required_argument, NULL, 't' },
		{ "seed", required_argument, NULL, 's' },
		{ "stdout", no_argument, NULL, 'o' },
		{ "json", no_argument, NULL, 'j' },
		{ "threads", required_argument, NULL, 'n' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, 0, NULL, 0 }
	};
	char *opts_short = ut_build_short_option_string(opts_long);

	int c, idx;
	char const *group_arg = NULL, *test_arg = NULL;
	while((c = getopt_long(argc, argv, opts_short, opts_long, &idx)) != -1) {
		switch(c) {
			case 'g': group_arg = optarg; break;
			case 't': test_arg = optarg; break;
			case 's': srand((unsigned int)atol(optarg)); break;
			case 'j': params->printer = ut_json_printer; break;
			case 'o': params->fp = stdout; break;
			case 'n': params->threads = atoi(optarg); break;
			case 'h': ut_print_help(); return(1);
			default: break;
		}
	}

	if(group_arg != NULL) {
		ut_modify_test_config_mark(group_arg, (void *)sorted_config, file_cnt);
	} else {
		ut_modify_test_config_all((void *)sorted_config, file_cnt);
	}

	if(test_arg != NULL) {
		ut_modify_test_config_mark(test_arg, (void *)sorted_test, test_cnt);
	} else {
		ut_modify_test_config_all((void *)sorted_test, test_cnt);
	}

	free(opts_short);
	return(0);
}

/**
 * @fn ut_propagate_config
 */
static inline
void ut_propagate_config(
	struct ut_s *test,
	int64_t test_cnt,
	struct ut_group_config_s *compd_config,
	uint64_t *sorted_file_idx,
	int64_t file_cnt)
{
	/* set index */
	for(uint64_t i = 0; i < file_cnt; i++) {
		for(uint64_t j = sorted_file_idx[i]; j < sorted_file_idx[i + 1]; j++) {
			test[j].index = i;
			test[j].succ = 0;		/* clear counters */
			test[j].fail = 0;
		}
	}
	return;
}

/**
 * @fn ut_run_test
 */
static
void ut_run_test(
	struct ut_s *test,
	struct ut_global_config_s const *gconf,
	struct ut_group_config_s const *compd_config)
{
	if(test->exec == 0) { return; }
	uint64_t index = test->index;

	/* initialize group context */
	void *gctx = NULL;
	if(compd_config[index].init != NULL && compd_config[index].clean != NULL) {
		gctx = compd_config[index].init(compd_config[index].params);
	}

	/* initialize local context */
	void *ctx = NULL;
	if(test->init != NULL && test->clean != NULL) {
		ctx = test->init(test->params);
	}

	/* run a test */
	test->fn(ctx, gctx, test, gconf, &compd_config[index]);

	/* cleanup contexts */
	if(test->init != NULL && test->clean != NULL) {
		test->clean(ctx);
	}
	if(compd_config[index].init != NULL && compd_config[index].clean != NULL) {
		compd_config[index].clean(gctx);
	}
	return;
}

/**
 * @fn ut_main_impl
 */
static
int ut_main_impl(int argc, char *argv[])
{
	/* dump symbol table */
	struct ut_nm_result_s *nm = ut_nm(argv[0]);

	if(nm == NULL) {
		return(1);
	}

	/* dump tests and configs */
	struct ut_s *test = ut_get_unittest(nm);
	struct ut_group_config_s *config = ut_get_ut_config(nm);

	/* sort by group, tag, line */
	ut_sort(test, config);
	struct ut_group_config_s *compd_config = ut_compensate_config(test, config);

	uint64_t test_cnt = ut_get_total_test_count(test);
	uint64_t *file_idx = ut_build_file_index(test);
	uint64_t file_cnt = ut_get_total_file_count(test);

	// printf("%" PRId64 "\n", file_cnt);

	/* topological sort by tag */
	for(uint64_t i = 0; i < file_cnt; i++) {
		// printf("toposort %" PRId64 ", %" PRId64 "\n", file_idx[i], file_idx[i + 1]);
		if(ut_toposort_by_tag(&test[file_idx[i]], file_idx[i + 1] - file_idx[i]) != 0) {
			fprintf(stderr, ut_color(UT_RED, "ERROR") ": failed to order tests. check if the depends_on options are sane.\n");
			return(1);
		}
	}

	// ut_dump_test(test);

	/* topological sort by group */
	if(ut_toposort_by_group(test, test_cnt, compd_config, file_idx, file_cnt) != 0) {
		fprintf(stderr, ut_color(UT_RED, "ERROR") ": failed to order tests. check if the depends_on options are sane.\n");
		return(1);
	}

	/* rebuild file idx */
	uint64_t *sorted_file_idx = ut_build_file_index(test);

	/* init default params */
	struct ut_global_config_s gconf = {
		.fp = stderr,
		.printer = ut_default_printer
	};

	/* modify config */
	if(ut_modify_test_config(argc, argv, &gconf, test, test_cnt, compd_config, file_cnt)) {
		return 1;
	}

	/* copy exec flag */
	ut_propagate_config(test, test_cnt, compd_config, sorted_file_idx, file_cnt);

	/* run tests */
	#ifdef _OPENMP
	omp_set_num_threads(gconf.threads);
	#pragma omp parallel for
	#endif
	for(uint64_t i = 0; i < test_cnt; i++) {
		ut_run_test(&test[i], &gconf, compd_config);
	}

	/* collect results */
	struct ut_result_s *res = calloc(sizeof(struct ut_result_s), file_cnt);
	for(uint64_t i = 0; i < test_cnt; i++) {
		uint64_t index = test[i].index;
		res[index].cnt++;
		res[index].succ += test[i].succ;
		res[index].fail += test[i].fail;
	}

	/* print results */
	gconf.printer.result(&gconf, compd_config, res, file_cnt);

	free(res);
	free(sorted_file_idx);
	free(file_idx);
	free(compd_config);
	free(test);
	free(config);
	free(nm);
	return(0);
}

static
int unittest_main(int argc, char *argv[])
{
	/* disable unittests if UNITTEST == 0 */
	#if UNITTEST != 0
		return(ut_main_impl(argc, argv));
	#else
		return(0);
	#endif
}

#if UNITTEST_ALIAS_MAIN != 0
int main(int argc, char *argv[])
{
	return(unittest_main(argc, argv));
}
#endif

/**
 * end of unittest.h
 */
