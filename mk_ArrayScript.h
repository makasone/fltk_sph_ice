#ifndef _MK_ARRAYSCRIPT_
#define _MK_ARRAYSCRIPT_

#include <iostream>

using namespace std;

class mk_ArrayScript
{
	//�z��̃T�C�Y��Ԃ�
public:
	template<typename TYPE, size_t SIZE> 
	static size_t array_length(const TYPE (&)[SIZE])
	{   
	    return SIZE;
	}

};

#endif;