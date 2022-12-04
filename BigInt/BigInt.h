#pragma once
#include <cstddef>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
using namespace std;

struct Big_Int
{
	size_t SIZE;
	unsigned long long* number;
	bool negative = false;

	Big_Int(size_t, unsigned long long); 
	Big_Int(string);
	Big_Int(unsigned long long);
	~Big_Int();
	Big_Int(const Big_Int&);

	Big_Int operator=(const Big_Int&);
	bool operator<(const Big_Int&) const;
	bool operator==(const Big_Int&) const;
	bool operator<=(const Big_Int&) const;
	friend ostream& operator<<(ostream& os, const Big_Int&);

	Big_Int operator+(const Big_Int&) const;
	Big_Int operator-(const Big_Int&) const;
	Big_Int operator*(const Big_Int&); // smart multiplication - chooses between Karatsuba & naive on itself;
	static Big_Int naive_multiplication(const Big_Int&,const Big_Int&);
	static Big_Int Karatsuba_multiplication(Big_Int, Big_Int);
	static const Big_Int _pow(const Big_Int&, const Big_Int&);
	Big_Int operator/(const Big_Int&);

	//Helper functions; unintended for manual usage
	void _squeeze();
	static void _to_equal_size(Big_Int&, Big_Int&);
	void _increaseByModulo();
	
	struct nested
	{
		size_t BASE = 10;
		size_t BASE_POWER = 3; // can't be more than 9
		size_t MODULO = pow(BASE,BASE_POWER);
		int Karatsuba_Power = 256;
	};
	inline static const nested insides; // C++17 or OLDER;


};