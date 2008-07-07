#ifndef CFX_FORTRAN_DEFS_H
#define CFX_FORTRAN_DEFS_H

//#include <g2c.h>

// external Fortran functions

#if defined (__cplusplus)
extern "C" {
#endif

    void openfile_(int *nerr,
                   char *filename,
                   int *whatfile,
                   int len_filename);

    void closefile_(int *nerr,
                    int *whatfile);

    void redsht_(int *dattyp, int *length, int *nerr,
                 char *what, char *where, int *when,
                 int *nsize,
                 int *iopt, int *ioptar,
                 float *rarr, int *iarr, char* carr,
                 int *larr, double *darr, char *sarr,
                 int len_what, int len_where, int len_carr);
    
    void readlong_(int *dattyp, int *nerr, char *what, char *where, 
                   int *when, int *nsize, int *iopt,
                   float *rarr, int *iarr, char* carr,
                   int *larr, double *darr, char *sarr,
                   int len_what, int len_where, int len_carr);

    void stri2s_(char *string, int* n, int len_string);

#if defined (__cplusplus)
}
#endif

/*
#include <g2c.h>

// external Fortran functions

#if defined (__cplusplus)
extern "C" {
#endif

    void openfile_(int *nerr,
                   char *filename,
                   int *whatfile,
                   ftnlen len_filename);

    void closefile_(int *nerr,
                    int *whatfile);

    void redsht_(int *dattyp, int *length, int *nerr,
                 char *what, char *where, int *when,
                 int *nsize,
                 int *iopt, int *ioptar,
                 float *rarr, int *iarr, char* carr,
                 int *larr, double *darr, char *sarr,
                 ftnlen len_what, ftnlen len_where, ftnlen len_carr);
    
    void readlong_(int *dattyp, int *nerr, char *what, char *where, 
                   int *when, int *nsize, int *iopt,
                   float *rarr, int *iarr, char* carr,
                   int *larr, double *darr, char *sarr,
                   ftnlen len_what, ftnlen len_where, ftnlen len_carr);

    void stri2s_(char *string, int* n, int len_string);

#if defined (__cplusplus)
}
#endif
*/

// Constants from cfd_io.h

//----------------------------------------------------------------- 
//  
//------ file units
//
//-----------------------------------------------------------------
//
#define   __unit_trm__ 1
#define   __unit_out__ 2
#define   __unit_err__ 3
#define   __unit_inp__ 4
//----------------------------------------------------------------- 
//  
//------ file names
//
//-----------------------------------------------------------------
//
#define   __def_file_name__ 'def'
#define   __res_file_name__ 'res'

#define   __out_file_name__ 'out'
#define   __out_file_number__ 50

#define   __job_file_name__ 'job'

#define   __hst_file_name__ 'hst'
#define   __phst_file_name__ 'phst'

#define   __stp_file_name__ 'stp'

#define   __msg_file_name__ 'english_msg.txt'

#define   __trn_file_name__ 'trn'

#define   __cgns_file_name__ 'cgns'

#define   __bnd_file_name__ 'bnd'

#define   __csv_file_name__ 'csv'

#define   __par_file_name__ 'par'

#define   __bak_file_name__ 'bak'

#define   __crash_file_name__ 'crash'

#define   __reread_file_name__ 'reread'

//-----------------------------------------------------------------
//
//------ file types
//
//-----------------------------------------------------------------

#define   __res_file_type__  0
#define   __trn_file_type__  1
#define   __bak_file_type__  2
#define   __cgns_file_type__ 3
#define   __bnd_file_type__  4
#define   __csv_file_type__  5

//-----------------------------------------------------------------
//
//------ outfile control level 
//
//-----------------------------------------------------------------

#define __brief__             0
#define __residual_only__     1
#define __verbose__           2


// Constants from cfx_fortran_io.h

//-------- define data type
//
#define            __real_data_type__ 1 
#define             __int_data_type__ 2 
#define            __char_data_type__ 3 
#define         __logical_data_type__ 4 
#define          __double_data_type__ 5 
#define          __string_data_type__ 6

#define          __io_offset_type__ INTEGER*8
#define          __io_string_type__ CHARACTER*(IO_BUFSIZE)

//
//-------- io error return code
//
#define                           __io_ok__  0
#define                     __io_open_err__  1 
#define                __io_file_not_open__  2
#define               __io_file_not_found__  3 
#define                    __io_disk_full__  4
#define                  __io_end_of_file__  5
#define            __io_permission_denied__  6
#define                    __io_close_err__  7
#define                     __io_seek_err__  8
#define                     __io_tell_err__  9
#define                     __io_read_err__ 10
#define                    __io_write_err__ 11  
#define                    __io_parse_err__ 12
#define              __io_unknown_dataset__ 13
#define                   __io_memory_err__ 14
#define                 __io_file_fmt_err__ 15
#define                __io_found_more_ds__ 16
#define                 __io_ds_not_found__ 17
#define                      __io_str_err__ 18
#define                     __io_mesg_err__ 19
#define                     __io_bad_file__ 20
#define                       __io_locked__ 21
#define                  __io_unknown_req__ 22 
#define                 __io_compress_err__ 23 
#define               __io_file_not_exist__ 24 
#define          __io_attribute_not_found__ 25
#define                          __io_err__ 26 
#define                  __io_convert_err__ 27
#define                 __io_memalloc_err__ 28

//
//-------- IOCNTF request list
//
#define             __io_open_primaryfile__  0
#define            __io_close_primaryfile__  1
#define           __io_open_secondaryfile__  2
#define          __io_close_secondaryfile__  3
#define                   __io_set_format__  4
#define               __io_delete_dataset__  5
#define                __io_write_dataset__  6
#define               __io_locate_dataset__  7
#define       __io_locate_primary_dataset__  8
#define       __io_delete_primary_dataset__  9
#define         __io_read_located_dataset__ 10
#define                  __io_read_header__ 11
#define                 __io_write_header__ 12
#define                    __io_read_data__ 13
#define                   __io_write_data__ 14
#define                __io_get_ds_number__ 15
#define                 __io_copy_dataset__ 16
#define             __io_copy_allvalid_ds__ 17
#define      __io_read_located_primary_ds__ 18
#define                 __io_write_endtag__ 19
#define                   __io_get_offset__ 20 
#define               __io_write_dslength__ 21
#define                  __io_goto_offset__ 22
#define           __io_set_index_location__ 23
#define                 __io_insert_index__ 24
#define           __io_get_index_location__ 25
#define                 __io_read_datatag__ 26
#define                __io_write_datatag__ 27
#define            __io_goto_primary_file__ 28
#define          __io_goto_secondary_file__ 29
#define                 __io_open_txtfile__ 30
#define                __io_close_txtfile__ 31
#define       __io_read_stri_from_txtfile__ 32
#define        __io_write_stri_to_txtfile__ 33
#define       __io_read_intr_from_txtfile__ 34
#define        __io_write_intr_to_txtfile__ 35
#define       __io_read_real_from_txtfile__ 36
#define        __io_write_real_to_txtfile__ 37
#define             __io_check_file_state__ 38
#define                __io_copy_when0_ds__ 39
#define          __io_insert_index_nosort__ 40 
#define                 __io_update_index__ 41
#define               __io_recreate_index__ 42
#define      __io_read_located_dataset_as__ 43
#define   __io_read_located_primary_ds_as__ 44
#define        __io_copy_ds_compress_when__ 45
#define           __io_set_default_format__ 46
//
//-------- action
//
#define                   __stop_if_failed__ 1
#define               __not_stop_if_failed__ 2
#define                   __warn_if_failed__ 3
#define                   __skip_if_failed__ 4

#endif // CFX_FORTRAN_DEFS_H
