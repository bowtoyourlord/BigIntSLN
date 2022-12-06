#include "BigInt.h"

const Big_Int MIN_FOR_KARATSUBA = Big_Int::_pow(10, Big_Int::insides.Karatsuba_Power);

Big_Int::Big_Int(size_t SIZE, unsigned long long fill_with)
{
	this->SIZE = SIZE;
	number = new unsigned long long[this->SIZE];
	for (size_t i = 0; i < this->SIZE; ++i)
		number[i] = fill_with;
}

Big_Int::Big_Int(string str)
{
	if (str[0] == '-')
	{
		negative = true;
		str.erase(0, 1);
	};
	if (str.size()< insides.BASE_POWER)
	{
		SIZE = 1;
		number = new unsigned long long[SIZE];
		number[0] = stoi(str);
	}
	else
	{
			size_t position = str.size();
			SIZE = str.size() / insides.BASE_POWER;
			if(str.size() % insides.BASE_POWER != 0) SIZE++;
			number = new unsigned long long[SIZE];
			string substring;
			for (size_t g = 0; g < SIZE-1; ++g)
			{
				position -= insides.BASE_POWER;
				substring = str.substr(position, insides.BASE_POWER);
				number[g] = stoi(substring);
			}
			substring = str.substr(0, position);
			number[SIZE - 1] = stoi(substring);
	}
}
Big_Int::Big_Int(unsigned long long other)
{
	SIZE = 1;
	if (other == 0)
	{
		number = new unsigned long long[SIZE];
		number[0] = 0;
	}
	else
	{
		size_t digits = log10(other) + 1;
		while (digits > insides.BASE_POWER)
		{
			digits -= insides.BASE_POWER;
			SIZE++;
		}
		number = new unsigned long long[SIZE];
		size_t i = 0;
		for (i; i < SIZE - 1; ++i)
		{
			number[i] = other % size_t(pow(insides.MODULO, (i + 1))); //cmath pow returns double ; explicit int cast required
			other /= size_t(pow(insides.MODULO, (i + 1)));
		}
		number[i] = other;
	};
};

Big_Int::~Big_Int()
{
	delete[] number;
}

Big_Int::Big_Int(const Big_Int& other)
{
	this->SIZE = other.SIZE;
	this->negative = other.negative;
	this->number = new unsigned long long[SIZE];
	for (auto i = 0; i < this->SIZE; i++)
	{
		this->number[i] = other.number[i];
	};
}

Big_Int Big_Int::operator=(const Big_Int& other)
{
	delete[] this->number;

	this->SIZE = other.SIZE;
	this->negative = other.negative;
	this->number = new unsigned long long[SIZE];
	for (auto i = 0; i < this->SIZE; i++)
	{
		this->number[i] = other.number[i];
	};

	return *this;
}

Big_Int Big_Int::operator+(const Big_Int& other) const
{
	if (this->negative && !(other.negative))
	{
		Big_Int tmp = *this;
		tmp.negative = false;
		Big_Int result = other - tmp;
		result._squeeze();
		return result;
	}
	else if (!(this->negative) && other.negative)
	{
		Big_Int tmp = other;
		tmp.negative = false;
		Big_Int result = *this - tmp;
		result._squeeze();
		return result;
	}

	Big_Int result(max(this->SIZE, other.SIZE) + 1, 0);
	if (this -> negative && other.negative) result.negative=true;

	size_t carry = 0;

	size_t i = 0;

	for (i; i < min(this->SIZE, other.SIZE); ++i)
	{
		size_t tmp = this->number[i] + other.number[i] + carry;
		carry = tmp / insides.MODULO;
		tmp %= insides.MODULO;
		result.number[i] = tmp;
	}
	if (this->SIZE == other.SIZE) result.number[i] = carry;
	else if (i == this->SIZE) for (i; i < other.SIZE;++i)
	{
		result.number[i] = other.number[i] + carry;
		carry = 0;
	}
	else for (i; i < this->SIZE;++i)
	{
		result.number[i] = this->number[i] + carry;
		carry = 0;
	}
	result._squeeze();
	return result;
}

Big_Int Big_Int::operator-(const Big_Int& other) const
{
	if ((* this) == other) return 0;
	if (other.negative)
	{
		Big_Int tmp = other;
		tmp.negative = false;
		Big_Int result = *this + tmp;
		result._squeeze();
		return result;
	}
	if (this->negative)
	{
		Big_Int tmp = *this;
		tmp.negative = false;
		tmp = tmp + other;
		tmp.negative = true;
		tmp._squeeze();
		return tmp;
	}
	if ((* this) < other)
	{
		Big_Int tmp = other - (*this);
		tmp.negative = true;
		tmp._squeeze();
		return tmp;
	}
	Big_Int result((this->SIZE) +1, 0);

	size_t carry = 0;

	size_t i = 0;

	for (i; i < other.SIZE; ++i)
	{
		size_t tmp = this->number[i] - other.number[i] - carry + insides.MODULO;
		if (tmp >= insides.MODULO)
		{
			tmp -= insides.MODULO;
			carry = 0;
		}
		else
		{
			carry = 1;
		}
		result.number[i] = tmp;
	}
	if (i < this->SIZE)
	{
		for (i; i < this->SIZE; ++i)
		{
			result.number[i] = this->number[i] - carry;
			carry = 0;
		}
	}
	result._squeeze();
	return result;


}

Big_Int Big_Int::operator*(const Big_Int& other)
{

	if (((*this) < MIN_FOR_KARATSUBA) || (other < MIN_FOR_KARATSUBA))
	{
		return naive_multiplication(*this, other);
	}
	else 
		return Karatsuba_multiplication(*this, other);

}
Big_Int Big_Int::naive_multiplication(const Big_Int& first,const Big_Int& other)
{
	if (first == 0 || other == 0)
		return 0;
	Big_Int result(first.SIZE + other.SIZE, 0);
	if (first.negative ^ other.negative) result.negative = true;
	{
		for (size_t i = 0; i < other.SIZE; ++i)
		{
			size_t carry = 0;
			for (size_t j = 0; j < first.SIZE; ++j)
			{
				result.number[i + j] += carry + first.number[j] * other.number[i];
				carry = result.number[i + j] / insides.MODULO;
				result.number[i + j] %= insides.MODULO;
			}
			result.number[i + first.SIZE] += carry;
		}
		result._squeeze();
		return result;
	}
}
Big_Int Big_Int::Karatsuba_multiplication(Big_Int first, Big_Int second)
{
	Big_Int result(first.SIZE + second.SIZE, 0);
	if (first.negative ^ second.negative) result.negative = true;
	_to_equal_size(first, second);
	size_t k = first.SIZE / 2;
	Big_Int x = _pow(10, insides.BASE_POWER*k);
	Big_Int x2 = _pow(x,2);

	Big_Int b(k, 0);
	for (size_t i = 0; i < k; ++i)
		b.number[i] = first.number[i];

	Big_Int a(first.SIZE - k, 0);
	for (size_t i = 0; i < a.SIZE; ++i)
		a.number[i] = first.number[i + k];

	Big_Int d(k, 0);
	for (size_t i = 0; i < k; ++i)
		d.number[i] = second.number[i];

	Big_Int c(second.SIZE - k, 0);
	for (size_t i = 0; i < c.SIZE; ++i)
		c.number[i] = second.number[i + k];

	Big_Int a_plus_b = a + b;
	Big_Int c_plus_d = c + d;
	Big_Int a_c = a * c;
	Big_Int b_d = b * d;
	Big_Int middle = a_plus_b * c_plus_d - a_c - b_d;
	result = (naive_multiplication(a_c, x2)) + naive_multiplication(middle, x) + b_d;
	result._squeeze();
	return result;
}
const Big_Int Big_Int::_pow(const Big_Int& base, const Big_Int& exponent)
{
	return exponent == 0 ? 1 : naive_multiplication(base,_pow(base, exponent - 1));
}
Big_Int Big_Int::operator/(const Big_Int& other)
{
	if (other == 0) return 0; // error
	if (*this == other) return 1;
	if (other == 1) return *this;
	Big_Int result(this->SIZE, 0); 
	if (this->negative && other.negative)
	{
		Big_Int tmp1 = *this;
		Big_Int tmp2 = other;
		tmp1.negative = false;
		tmp2.negative = false;
		result = tmp1 / tmp2;
		result._squeeze();
		return result;
	}
	if (this->negative) 
	{
		Big_Int tmp = *this;
		tmp.negative = false;
		result = (tmp / other);
		result.negative = true;
		result._squeeze();
		return result;
	}
	if (other.negative)
	{
		Big_Int tmp = other;
		tmp.negative = false;
		result = (*this / tmp);
		result.negative = true;
		result._squeeze();
		return result;
	}
	if (this->negative && other.negative) if (other < *this) return 0;
	else if (*this < other) return 0;

	Big_Int current("0");

	for (auto i = this->SIZE -1; i != -1; --i)
	{
		current.number[0] = this->number[i];
		unsigned long long x = 0, l = 0, r = insides.MODULO;
		while (l <= r)
		{
			auto m = (l + r) / 2;
			Big_Int minus = Big_Int(m) * other;
			if (minus <= current)
			{
				x = m;
				l = m + 1;
			}
			else
				r = m - 1;
		}
		result.number[i] = x;
		current = current - (Big_Int(x)*other);
		current._increaseByModulo();
	}
	result._squeeze();
	return result;

}

void Big_Int::_squeeze()
{
	if ((*this) == Big_Int("0")) return;
	size_t counter = 0;
	for (size_t i = this->SIZE - 1; i != -1; --i)
	{
		if (this->number[i] == 0) counter++;
		else break;
	}
	if (counter == 0) return;
	if (counter == SIZE)
	{
		*this = Big_Int("0");
		return;
	};
	unsigned long long* tmp = new unsigned long long[this->SIZE - counter];
	for (size_t i = 0; i < this->SIZE - counter; ++i) tmp[i] = number[i];
	delete[] this->number;
	SIZE -= counter;
	this->number = new unsigned long long[SIZE];
	for (size_t i = 0; i < this->SIZE; ++i) number[i] = tmp[i];
	delete[] tmp;
}

void Big_Int::_to_equal_size(Big_Int& first, Big_Int& other)
{
	if (first.SIZE == other.SIZE) return;
	Big_Int tmp(max(first.SIZE, other.SIZE), 0);
	if (first.SIZE < other.SIZE) for (size_t i = 0; i < first.SIZE; ++i)
	{
		tmp.number[i] = first.number[i];
		first = tmp;
	}
	else for (size_t i = 0; i < other.SIZE; ++i)
	{
		tmp.number[i] = other.number[i];
		other = tmp;
	}
}

void Big_Int::_increaseByModulo()
{

	unsigned long long* tmp = new unsigned long long[this->SIZE + 1];
	for (size_t i = this->SIZE - 1; i != -1; --i)
		tmp[i + 1] = this->number[i];
	delete[] this->number;
	tmp[0] = 0;
	this->number = new unsigned long long[++SIZE];
	for (size_t i = 0; i <this->SIZE; ++i)
		this->number[i] = tmp[i];
	delete[] tmp;
}


bool Big_Int::operator<(const Big_Int& other) const
{
	if (this->negative && other.negative)
	{
		if (this->SIZE < other.SIZE) return false;
		else if (other.SIZE < this->SIZE) return true;
		else
		{
			for (size_t i = SIZE - 1; i != -1; --i)
			{
				if (this->number[i] < other.number[i]) return false;
				else if (other.number[i] < this->number[i]) return true;
			}
			return false;
		}
	}
	else if (!(other.negative) && this->negative) return true;
	else if (other.negative && !(this->negative)) return false;
	else
	{
		if (this->SIZE < other.SIZE) return true;
		else if (other.SIZE < this->SIZE) return false;
		else 
		{
			for (size_t i = SIZE - 1; i != -1; --i)
			{
				if (this->number[i] < other.number[i]) return true;
				else if (other.number[i] < this->number[i]) return false;
			}
			return false;
		}
	}
}

bool Big_Int::operator==(const Big_Int& other) const
{
	if (this->SIZE == other.SIZE && this->negative == other.negative)
	{
		for (size_t i = 0; i < SIZE; ++i)
		{
			if (this->number[i] != other.number[i]) return false;
		}
		return true;
	}
	else return false;

}

bool Big_Int::operator<=(const Big_Int& other) const
{
	return ((*this) < (other) || (*this) == (other));
}

ostream& operator<<(ostream& os, const Big_Int& obj)
{
	if (obj == 0) return os << "0";
	string result;
	if (obj.negative) result.push_back('-');
	for (size_t i = obj.SIZE - 1; i != -1; --i)
	{
		if (obj.number[i] == 0)
		{
			string tmp = "0";
			tmp.insert(0, Big_Int::insides.BASE_POWER-1, '0');
			result += tmp;
		}
		else if (i == obj.SIZE - 1)
		{
			string tmp = to_string(obj.number[i]);
			result += tmp;
		}
		else 
		{
			string tmp = to_string(obj.number[i]);
			if (tmp.size() < Big_Int::insides.BASE_POWER)
			{
				string nulls = "0";
				nulls.insert(0, Big_Int::insides.BASE_POWER - tmp.size() - 1, '0');
				tmp = nulls + tmp;
			}
				

			result += tmp;
		}
	}
	return os << result;

}
