#ifndef CASATOOLS_DATA_OPT_H_
#define CASATOOLS_DATA_OPT_H_
#include <stdexcept>

namespace casatools {

    template <typename T> class OptionValue {
    public:
        const T &get( ) {
            if ( ! has_value_ )
                throw std::runtime_error("attempt to retrieve value from none (option)");
            return value_;
        }
        bool has_value( ) { return has_value_; }
    private:
        template<typename U> friend class opt;
        OptionValue( ) : has_value_(false) { }
        OptionValue( const T &val ): has_value_(true), value_(val) { }
        bool has_value_;
        T value_;
    };

    template <typename T> class opt {
    public:
        static OptionValue<T> some(const T& val) { return OptionValue<T>(val); }
        static OptionValue<T> none( ) { return OptionValue<T>( ); }
    };

}

#endif
