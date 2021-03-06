/*
	Author: John Grime, The University of Chicago.
*/

#if !defined( UTIL_STRING_DEFINED )
#define UTIL_STRING_DEFINED

#include <cstdlib>
#include <string.h>

#include <string>
#include <vector>
#include <map>

namespace Util
{

//
// As a class, to avoid unused method warnings from static "free" routines in a a namespace.
//
class String
{
	protected:
	
		//
		// Check if character 'test' is in test_characters (assumes "test_characters" null-terminated!).
		//
		static bool is_in( char test, const char * test_characters )
		{
			if( test_characters == nullptr ) return false;

			while( *test_characters != '\0' )
			{
				if( test == *(test_characters++) ) return true;
			}
			return false;
		}

		//
		// Check if "check" starts with string "str"; used in variable expansion routines.
		//
		static bool starts_with( const char *str, const char *check, int check_len )
		{
			if( str[0] == '\0' || check[0] == '\0' ) return false;
	
			for( int i=0; i<check_len && str[i] != '\0'; i++ )
			{
				if( str[i] != check[i] ) return false;
			}
			return true;
		}

	public:

		enum class ReturnValue { ERROR, OK };

		//
		// Convert a character sequence into an integer or floating point type
		//
		template< typename T >
		static ReturnValue ToInteger( const char *str, T& result, int base = 10 )
		{
			if( str == nullptr || base < 2 ) return ReturnValue::ERROR;

			char *endptr;
			result = strtoll( str, &endptr, base );
			if( (endptr==str) || (*endptr!='\0') ) return ReturnValue::ERROR;
			return ReturnValue::OK;
		}
		template< typename T >
		static ReturnValue ToInteger( const std::string& str, T& result, int base = 10 )
		{
			return ToInteger( str.c_str(), result, base );
		}

		template< typename T >
		static ReturnValue ToReal( const char *str, T& result )
		{
			if( str == nullptr ) return ReturnValue::ERROR;

			char *endptr;
			result = strtod( str, &endptr );
			if( (endptr==str) || (*endptr!='\0') ) return ReturnValue::ERROR;
			return ReturnValue::OK;
		}
		template< typename T >
		static ReturnValue ToReal( const std::string& str, T& result )
		{
			return ToReal( str.c_str(), result );
		}

		template< typename T >
		static ReturnValue ToNumber( const char *str, T& result, int base = 10 )
		{
			if( std::is_integral<T>::value ) return ToInteger( str, result, base );
			else if( std::is_floating_point<T>::value ) return ToReal( str, result );
			return ReturnValue::ERROR;
		}
		template< typename T >
		static ReturnValue ToNumber( const std::string& str, T& result, int base = 10 )
		{
			return ToNumber( str.c_str(), result, base );
		}

		//
		// Extract a set of integers from an appropriately formatted string, using "element_sep"
		// and "range_sep" strings. E.g. : element_sep == "," and range_sep == "-" with the following
		// input:
		//
		//  12,13,14,20-22,90
		//
		// gives the list of integers: 12 13 14 20 21 22 90
		//
		// NOTES:
		// 1. Ranges are INCLUSIVE of the specified start and end value!
		// 2. "results" vector IS NOT CLEARED before new values appended, to allow progressive generation.
		//
		template <typename T>
		static ReturnValue ToIntegers( const char *str, std::vector<T> &results, const char *element_sep, const char *range_sep, int base = 10 )
		{
			std::vector< std::string > tokens, subtokens;
			int start = -1, stop = -1; // will ALWAYS be set, but stop compiler moaning

			// tokens are either individual numbers, or numbers separated by the range token.
			Tokenize( str, tokens, element_sep );
			for( const auto& token : tokens )
			{
				// If splitting this token on the range character gives two tokens, assume they are range start & stop
				if( Tokenize( token, subtokens, range_sep ) == 2 )
				{
					if( ToInteger( subtokens[0], start, base ) == ReturnValue::ERROR ) return ReturnValue::ERROR;
					if( ToInteger( subtokens[1], stop, base ) == ReturnValue::ERROR ) return ReturnValue::ERROR;
					if( stop < start ) return ReturnValue::ERROR;
				}
				else
				{
					if( ToInteger( token, start, base ) == ReturnValue::ERROR ) return ReturnValue::ERROR;
					stop = start;
				}
				
				for( int value=start; value<=stop; value++ ) results.push_back( value );
			}
			return ReturnValue::OK;
		}
		template <typename T>
		static ReturnValue ToIntegers( const std::string& str, std::vector<T> &results, const char *element_sep, const char *range_sep, int base = 10 )
		{
			return ToIntegers( str.c_str(), results, element_sep, range_sep, base );
		}

		//
		// This is slow, as we repeatedly append to an std::string, and then copy it into "results". This approach
		// does mean we don't need any nasty fixed-size intermediate buffers etc, so it's hopefully safer.
		//
		static int Tokenize( const char * source, std::vector< std::string > &results, const char *delimiters )
		{
			std::string temp;
	
			results.clear(); // BEFORE error return test ...

			if( source == nullptr || delimiters == nullptr ) return 0;

			size_t src_len = strlen( source );
	
			for( size_t i=0; i<src_len; i++ )
			{
				if( is_in( source[i], delimiters ) == true )
				{
					if( temp.size() > 0 ) results.push_back( temp );
					temp.clear();
				}
				else
				{
					temp += source[i];
				}
			}
			// If there's a trailing token, add to the results vector
			if( temp.size() > 0  ) results.push_back( temp );
	
			return (int)results.size();
		}
		static int Tokenize( const std::string& source, std::vector< std::string > &results, const char *delimiters )
		{
			return Tokenize( source.c_str(), results, delimiters );
		}
		
		//
		// Strip defined whitespace from start and end of string. MODIFIES "str"!
		//
		static ReturnValue Strip( char *str, const char *whitespace )
		{
			int str_length, start, end;

			if( str == nullptr || whitespace == nullptr ) return ReturnValue::OK;

			str_length = strlen( str );

			// get index of first non-whitespace character in source string
			start = -1;
			for( int i=0; i<str_length && start == -1; i++ )
			{
				if( is_in( str[i], whitespace ) == false ) start = i;
			}

			// get index of last non-whitespace character in source string
			end = -1;
			for( int i=str_length-1; i>=0 && end == -1; i-- )
			{
				if( is_in( str[i], whitespace ) == false ) end = i;
			}

			// note that we could have an empty string, so check for that.
			if( start == -1 || end == -1 || end < start )
			{
				str[0] = '\0';
			}
			else
			{
				int j = 0;
				for( int i=start; i<=end; i++ ) str[j++] = str[i];
				str[j] = '\0';
			}

			return ReturnValue::OK;
		}
		// No std::string version of the above: I'm not sure how implementtions would
		// handle us potentially inserting a null char into a std::string? Would this
		// e.g. invaidate the size() method of std::string?

		//
		// Convert arbitrary POD into a character bitstring. Separates groups of "unit_size_bits"
		// with a space for easier reading. WARNING: assumes character buffer is large enough!
		// Here we assume a byte is 8 bits.
		//
		template <typename T>
		static void FromBinary( T data, char *cbuf, int unit_size_bits = 8 )
		{
			T one = 1;
			int nbits = (int)sizeof(T)*8, upto = 0;
	
			if( cbuf == nullptr ) return;

			for( int i=nbits-1; i>=0; i-- )
			{
				cbuf[upto++] = ( data & (one<<i) ) ? '1' : '0';
				if( i > 0 && i%unit_size_bits == 0 ) cbuf[upto++] = ' ';
			}
			cbuf[upto] = '\0';
		}

		//
		// Given an input string "in_str" (and a map of variables "vars"), create an output string "out_str" which
		// has the variables expanded into their values as per the "vars" map.
		//
		// lbracket: "open tag" for variable name
		// rbracket: "close tag" for variable name
		//
		// eg : lbracket == "${", rbracket == "}" would detect variables such as "${name}" and replace "${name}"
		// with the string found under the key "name" in "vars".
		//
		static ReturnValue ExpandVariables( const char *in_str, std::string &out_str, const char *lbracket, const char *rbracket, std::map<std::string,std::string> &vars )
		{
			int lb_len, rb_len, state;
			const char *in_str_upto;
			std::string var_name;
	
			if( in_str == nullptr || lbracket == nullptr || rbracket == nullptr ) return ReturnValue::ERROR;
	
			lb_len = strlen( lbracket );
			rb_len = strlen( rbracket );
	
			var_name = "";
			state = 0; // 0 == not in variable name parse, 1 == in variable name parse
			in_str_upto = in_str;
	
			out_str = "";
	
			while( *in_str_upto != '\0' )
			{
				//
				// Test to see if the current character starts a special variable bracketing string ...
				//
				bool hit_lbracket = starts_with( in_str_upto, lbracket, lb_len );
				bool hit_rbracket = starts_with( in_str_upto, rbracket, rb_len );
		
				//
				// Currently in a variable?
				//
				// - if hit rbracket: look up variable and replace in output string
				// - if hit lbracket: error - can't nest variable names at the moment
				// - otherwise: copy character into variable name
				//
				if( state == 1 )
				{
					if( hit_rbracket )
					{
						// Empty variable name?
						if( var_name == "" ) return ReturnValue::ERROR;
				
						// Variable not found?
						auto it = vars.find( var_name );
						if( it == end(vars) ) return ReturnValue::ERROR;
				
						out_str += it->second;
				
						state = 0;
						var_name = "";
						in_str_upto += rb_len;
					}
					else if( hit_lbracket )
					{
						// Unexpected variable start!
						return ReturnValue::ERROR;
					}
					else
					{
						var_name += *(in_str_upto++);
					}
				}
				//
				// Not in variable; just copy characters, unless we hit a variable bracketing string.
				//
				else
				{
					if( hit_lbracket )
					{
						state = 1;
						var_name = "";
						in_str_upto += lb_len;
					}
					else
					{
						out_str += *(in_str_upto++);
					}
				}
			}

			// String ended with open variable?
			if( state != 0 ) return ReturnValue::ERROR;
	
			return ReturnValue::OK;
		}
		static ReturnValue ExpandVariables( const std::string& in_str, std::string &out_str, const char *lbracket, const char *rbracket, std::map<std::string,std::string> &vars )
		{
			return ExpandVariables( in_str.c_str(), out_str, lbracket, rbracket, vars );
		}

		static const char* ReturnValueAsString( ReturnValue r )
		{
			switch( r )
			{
				case ReturnValue::ERROR:
					return "ERROR";
				break;

				case ReturnValue::OK:
					return "OK";
				break;
			}
			return "UNKNOWN";
		}
};

}

#endif
