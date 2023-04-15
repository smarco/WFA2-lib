/*
 * Starting from C++17, the boost filesystem library is supported.
 * However, gcc versions 7 and older only implemented the functionality as
 * experimental. 
 */

#ifdef __GNUC__
    #if __GNUC__ < 8
        #define FILESYSTEM_EXPERIMENTAL
    #endif
#endif

#ifdef FILESYSTEM_EXPERIMENTAL
    #include <experimental/filesystem>
    using namespace std::experimental;
#else
    #include <filesystem>
#endif
