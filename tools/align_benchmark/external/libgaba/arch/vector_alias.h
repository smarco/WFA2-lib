
/**
 * @file vector_alias.h
 *
 * @brief make alias to vector macros
 */
#if !defined(_NVEC_ALIAS_PREFIX) || !defined(_WVEC_ALIAS_PREFIX)
#error "_NVEC_ALIAS_PREFIX and _WVEC_ALIAS_PREFIX must be defined when alias.h is included."
#else /* _VECTOR_ALIAS_PREFIX */

#define nvec_prefix 					_NVEC_ALIAS_PREFIX			/* prefix for narrow variables */
#define wvec_prefix 					_WVEC_ALIAS_PREFIX			/* prefix for wide variables */

/* join macros */
#define _vec_alias_join2_intl(a,b)		a##b
#define _vec_alias_join2(a,b)			_vec_alias_join2_intl(a,b)

#define _vec_alias_join3_intl(a,b,c)	a##b##_##c
#define _vec_alias_join3(a,b,c)			_vec_alias_join3_intl(a,b,c)


/* types */
#define nvec_t 			_vec_alias_join2(nvec_prefix, _t)
#define nvec_mask_t		_vec_alias_join2(nvec_prefix, _mask_t)
#define nvec_masku_t	_vec_alias_join2(nvec_prefix, _masku_t)

#define wvec_t 			_vec_alias_join2(wvec_prefix, _t)
#define wvec_mask_t		_vec_alias_join2(wvec_prefix, _mask_t)
#define wvec_masku_t	_vec_alias_join2(wvec_prefix, _masku_t)

/* address */
#define _pv_n 			_vec_alias_join2(_pv_, nvec_prefix)
#define _pv_w 			_vec_alias_join2(_pv_, wvec_prefix)

/* load and store */
#define _load_n			_vec_alias_join2(_load_, nvec_prefix)
#define _loadu_n		_vec_alias_join2(_loadu_, nvec_prefix)
#define _store_n		_vec_alias_join2(_store_, nvec_prefix)
#define _storeu_n		_vec_alias_join2(_storeu_, nvec_prefix)

#define _load_w			_vec_alias_join2(_load_, wvec_prefix)
#define _loadu_w		_vec_alias_join2(_loadu_, wvec_prefix)
#define _store_w		_vec_alias_join2(_store_, wvec_prefix)
#define _storeu_w		_vec_alias_join2(_storeu_, wvec_prefix)

/* broadcast */
#define _set_n			_vec_alias_join2(_set_, nvec_prefix)
#define _zero_n			_vec_alias_join2(_zero_, nvec_prefix)
#define _seta_n			_vec_alias_join2(_seta_, nvec_prefix)
#define _swap_n			_vec_alias_join2(_swap_, nvec_prefix)

#define _set_w			_vec_alias_join2(_set_, wvec_prefix)
#define _zero_w			_vec_alias_join2(_zero_, wvec_prefix)
#define _seta_w			_vec_alias_join2(_seta_, wvec_prefix)
#define _swap_w			_vec_alias_join2(_swap_, wvec_prefix)

/* logics */
#define _not_n			_vec_alias_join2(_not_, nvec_prefix)
#define _and_n			_vec_alias_join2(_and_, nvec_prefix)
#define _or_n			_vec_alias_join2(_or_, nvec_prefix)
#define _xor_n			_vec_alias_join2(_xor_, nvec_prefix)
#define _andn_n			_vec_alias_join2(_andn_, nvec_prefix)

#define _not_w			_vec_alias_join2(_not_, wvec_prefix)
#define _and_w			_vec_alias_join2(_and_, wvec_prefix)
#define _or_w			_vec_alias_join2(_or_, wvec_prefix)
#define _xor_w			_vec_alias_join2(_xor_, wvec_prefix)
#define _andn_w			_vec_alias_join2(_andn_, wvec_prefix)

/* arithmetics */
#define _add_n			_vec_alias_join2(_add_, nvec_prefix)
#define _sub_n			_vec_alias_join2(_sub_, nvec_prefix)
#define _adds_n			_vec_alias_join2(_adds_, nvec_prefix)
#define _subs_n			_vec_alias_join2(_subs_, nvec_prefix)
#define _max_n			_vec_alias_join2(_max_, nvec_prefix)
#define _min_n			_vec_alias_join2(_min_, nvec_prefix)

#define _add_w			_vec_alias_join2(_add_, wvec_prefix)
#define _sub_w			_vec_alias_join2(_sub_, wvec_prefix)
#define _adds_w			_vec_alias_join2(_adds_, wvec_prefix)
#define _subs_w			_vec_alias_join2(_subs_, wvec_prefix)
#define _max_w			_vec_alias_join2(_max_, wvec_prefix)
#define _min_w			_vec_alias_join2(_min_, wvec_prefix)

/* shuffle */
#define _shuf_n			_vec_alias_join2(_shuf_, nvec_prefix)
#define _shuf_w			_vec_alias_join2(_shuf_, wvec_prefix)

/* blend */
#define _sel_n			_vec_alias_join2(_sel_, nvec_prefix)
#define _sel_w			_vec_alias_join2(_sel_, wvec_prefix)

/* compare */
#define _eq_n			_vec_alias_join2(_eq_, nvec_prefix)
#define _lt_n			_vec_alias_join2(_lt_, nvec_prefix)
#define _gt_n			_vec_alias_join2(_gt_, nvec_prefix)

#define _eq_w			_vec_alias_join2(_eq_, wvec_prefix)
#define _lt_w			_vec_alias_join2(_lt_, wvec_prefix)
#define _gt_w			_vec_alias_join2(_gt_, wvec_prefix)

/* insert and extract */
#define _ins_n			_vec_alias_join2(_ins_, nvec_prefix)
#define _ext_n			_vec_alias_join2(_ext_, nvec_prefix)

#define _ins_w			_vec_alias_join2(_ins_, wvec_prefix)
#define _ext_w			_vec_alias_join2(_ext_, wvec_prefix)

/* shift */
#define _bsl_n			_vec_alias_join2(_bsl_, nvec_prefix)
#define _bsr_n			_vec_alias_join2(_bsr_, nvec_prefix)
#define _shl_n			_vec_alias_join2(_shl_, nvec_prefix)
#define _shr_n			_vec_alias_join2(_shr_, nvec_prefix)
#define _sal_n			_vec_alias_join2(_sal_, nvec_prefix)
#define _sar_n			_vec_alias_join2(_sar_, nvec_prefix)

#define _bsl_w			_vec_alias_join2(_bsl_, wvec_prefix)
#define _bsr_w			_vec_alias_join2(_bsr_, wvec_prefix)
#define _shl_w			_vec_alias_join2(_shl_, wvec_prefix)
#define _shr_w			_vec_alias_join2(_shr_, wvec_prefix)
#define _sal_w			_vec_alias_join2(_sal_, wvec_prefix)
#define _sar_w			_vec_alias_join2(_sar_, wvec_prefix)

/* mask */
#define _mask_n			_vec_alias_join2(_mask_, nvec_prefix)
#define _mask_w			_vec_alias_join2(_mask_, wvec_prefix)

/* horizontal max */
#define _hmax_n			_vec_alias_join2(_hmax_, nvec_prefix)
#define _hmax_w			_vec_alias_join2(_hmax_, wvec_prefix)

/* broadcast */
#define _from_v16i8_n	_vec_alias_join2(_from_v16i8_, nvec_prefix)
#define _from_v32i8_n	_vec_alias_join2(_from_v32i8_, nvec_prefix)
#define _to_v16i8_n		_vec_alias_join2(_to_v16i8_, nvec_prefix)
#define _to_v32i8_n		_vec_alias_join2(_to_v32i8_, nvec_prefix)

/* convert */
#define _cvt_n_w		_vec_alias_join3(_cvt_, nvec_prefix, wvec_prefix)
#define _cvt_w_n		_vec_alias_join3(_cvt_, wvec_prefix, nvec_prefix)

/* print */
#define _print_n		_vec_alias_join2(_print_, nvec_prefix)
#define _print_w		_vec_alias_join2(_print_, wvec_prefix)

#endif /* _VECTOR_ALIAS_PREFIX */
/**
 * end of vector_alias.h
 */
