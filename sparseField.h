#ifndef SPARSEFIELD_H
#define SPARSEFIELD_H

#include <map>
#include <itkSimpleFastMutexLock.h>


template< class TElement, class TIndex = int >
class sparsefield {
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename TIndex size_type;
	typedef typename TIndex difference_type;
	typedef typename TElement value_type;
	typedef typename const TElement const_reference;

	typedef std::map< TIndex, TElement > ContainerType;

	class iter_base {
	public:
		iter_base( ):m_Position(-1) { }
		iter_base( sparsefield &field, TIndex pos ): m_Field( &field), m_Position( pos ) {}
		iter_base( sparsefield &field ): m_Field( &field), m_Position( 0 ) {}
	protected:
		TIndex m_Position;
		sparsefield *m_Field;
	};

	class reference: protected iter_base {	// reference 
	public:
		reference()	{ }
		reference(const iter_base& _Right): iter_base(_Right) {	// construct with base
		}
		reference& operator=(const reference& _Right) {	// assign reference _Right to bit
			return (*this = TElement(_Right));
		}
		reference& operator=(TElement _Val)	{	// assign _Val to bit
			m_Field->Set( _Val, m_Position );
			return (*this);
		}
		operator TElement() const {	// test if bit is set
			return m_Field->Get( m_Position );
		}
	};

	typedef typename reference *pointer;


	class const_iterator: public iter_base {
	public:
		typedef sparsefield baseclass;
		typedef typename baseclass::iterator_category iterator_category;
		typedef typename baseclass::value_type value_type;
		typedef typename baseclass::difference_type difference_type;
		typedef typename baseclass::reference reference;
		typedef typename baseclass::const_reference const_reference;
		typedef typename baseclass::pointer pointer;
		const_iterator(): iter_base( ) {}
		const_iterator( sparsefield &field, TIndex pos ): iter_base( field, pos ) {}
		const_iterator( sparsefield &field ): iter_base( field) {}
		const_reference operator*() const {	
			return const_reference( (*m_Field)[ m_Position ] );
		}
/*
		_Ctptr operator->() const {
			return (&**this);
		}
*/
		const_iterator& operator++() {	// preincrement
			++m_Position;
			return (*this);
		}

		const_iterator operator++(int) {	// postincrement
			const_iterator _Tmp = *this;
			++*this;
			return (_Tmp);
		}

		const_iterator& operator--() {	// predecrement
			--m_Position;
			return (*this);
		}

		const_iterator operator--(int) {	// postdecrement
			const_iterator _Tmp = *this;
			--*this;
			return (_Tmp);
		}

		const_iterator& operator+=(difference_type _Off) {	// increment by integer
			m_Position += _Off;
			return (*this);
		}

		const_iterator operator+(difference_type _Off) const {	// return this + integer
			const_iterator _Tmp = *this;
			return (_Tmp += _Off);
		}

		const_iterator& operator-=(difference_type _Off) {	// decrement by integer
			return (*this += -_Off);
		}

		const_iterator operator-(difference_type _Off) const {	// return this - integer
			const_iterator _Tmp = *this;
			return (_Tmp -= _Off);
		}

		difference_type operator-(const const_iterator& _Right) const {	// return difference of iterators
			return (m_Position - _Right.m_Position);
		}

		const_reference operator[](difference_type _Off) const {	// subscript
			return (*(*this + _Off));
		}

		bool operator==(const const_iterator& _Right) const	{	// test for iterator equality
			return (m_Position == _Right.m_Position);
		}

		bool operator!=(const const_iterator& _Right) const	{	// test for iterator inequality
			return (!(*this == _Right));
		}

		bool operator<(const const_iterator& _Right) const	{	// test if this < _Right
			return (m_Position < _Right.m_Position);
		}

		bool operator>(const const_iterator& _Right) const	{	// test if this > _Right
			return (_Right < *this);
		}

		bool operator<=(const const_iterator& _Right) const	{	// test if this <= _Right
			return (!(_Right < *this));
		}

		bool operator>=(const const_iterator& _Right) const	{	// test if this >= _Right
			return (!(*this < _Right));
		}

		friend const_iterator operator+(difference_type _Off,
			const const_iterator& _Right) {	// return iterator + integer
			return (_Right + _Off);
		}
		protected:
		};

		// CLASS iterator
	class iterator;
	friend class iterator;

	class iterator: public const_iterator {	// iterator for mutable sparsefield
	public:
		typedef sparsefield baseclass;
		typedef typename baseclass::iterator_category iterator_category;
		typedef typename baseclass::value_type value_type;
		typedef typename baseclass::difference_type difference_type;
		typedef typename baseclass::reference reference;
		typedef typename baseclass::const_reference const_reference;
		typedef typename baseclass::pointer pointer;
		iterator(): const_iterator() {}
		iterator( sparsefield &field, TIndex pos ): const_iterator( field, pos ) {}
		iterator( sparsefield &field ): const_iterator( field)) {}
		reference operator*() const	{	// return designated object
			return reference( *this);
		}
/*
		_Tptr operator->() const
			{	// return pointer to class object
			return (&**this);
			}
*/
		iterator& operator++()	{	// preincrement
			++m_Position;
			return (*this);
		}

		iterator operator++(int) {	// postincrement
			iterator _Tmp = *this;
			++*this;
			return (_Tmp);
		}

		iterator& operator--()	{	// predecrement
			--m_Position;
			return (*this);
		}

		iterator operator--(int) {	// postdecrement
			iterator _Tmp = *this;
			--*this;
			return (_Tmp);
		}

		iterator& operator+=(difference_type _Off)	{	// increment by integer
			m_Position += _Off;
			return (*this);
		}

		iterator operator+(difference_type _Off) const {	// return this + integer
			iterator _Tmp = *this;
			return (_Tmp += _Off);
		}

		iterator& operator-=(difference_type _Off)	{	// decrement by integer
			return (*this += -_Off);
		}

		iterator operator-(difference_type _Off) const	{	// return this - integer
			iterator _Tmp = *this;
			return (_Tmp -= _Off);
		}

		difference_type operator-(const const_iterator& _Right) const {	// return difference of iterators
			return ((const_iterator)*this - _Right);
		}

		reference operator[](difference_type _Off) const {	// subscript
			return (*(*this + _Off));
		}

		friend iterator operator+(difference_type _Off,	const iterator& _Right) 	{	// return iterator + integer
			return (_Right + _Off);
		}
	};
	

	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;


	sparsefield(): m_Size(0) { }

	explicit sparsefield(size_type _Count)
		: m_Size( _Count ) {
		if (_Count)	m_Container[ 0 ] = TElement();
	}

	sparsefield(size_type _Count, const TElement& _Val)
		: m_Size( _Count ) {
			m_Container[ 0 ] = val;	
		}

	sparsefield(const sparsefield& _Right)
		: m_Size(_Right.m_Size), m_Container( _Right.m_Container)	{}

	template<class _Iter>
		sparsefield(_Iter _First, _Iter _Last)
		{	
			m_Size = 0;
			if (_First!=_Last) {
				TElement old = *(_First++);
				m_Container[m_Size++] = old;
				while(_First!=_Last) {
					TElement curr = *_First;
					if (curr != old) {
						m_Container[ m_Size ] = curr;
						old = curr;
					}
					++_First;
					++m_Size;
				}
			}
		}
/*
	template<class _Iter>
		void _Construct(_Iter _Count, _Iter _Val, _Int_iterator_tag)
		{	// initialize with _Count * _Val
		size_type _Size = (size_type)_Count;
		_Construct_n(_Size, _Val);
		}
*/
	sparsefield& operator=(const sparsefield& _Right) {	// assign _Right
		if (this != &_Right) {	// worth doing
			m_Size = _Right.m_Size;
			m_Container = _Right.m_Container;
		}
		return (*this);
	}

	iterator begin()
		{	// return iterator for beginning of mutable sequence
		return (iterator(*this,0));
		}

	const_iterator begin() const
		{	// return iterator for beginning of nonmutable sequence
		return (const_iterator(*this,0));
		}

	iterator end()
		{	// return iterator for end of mutable sequence
		return (iterator(*this,m_Size));
		}

	const_iterator end() const
		{	// return iterator for end of nonmutable sequence
		return (const_iterator(*this,m_Size));
		}

	reverse_iterator rbegin()
		{	// return iterator for beginning of reversed mutable sequence
		return (reverse_iterator(end()));
		}

	const_reverse_iterator rbegin() const
		{	// return iterator for beginning of reversed nonmutable sequence
		return (const_reverse_iterator(end()));
		}

	reverse_iterator rend()
		{	// return iterator for end of reversed mutable sequence
		return (reverse_iterator(begin()));
		}

	const_reverse_iterator rend() const
		{	// return iterator for end of reversed nonmutable sequence
		return (const_reverse_iterator(begin()));
		}

	void resize(size_type _Newsize)	{	// determine new length, padding with _Ty() elements as needed
		m_MutexLock.Lock();
		if (!m_Size) m_Container[0] = TElement();
		else 	m_Container.erase( m_Container.lower_bound( _Newsize ), m_Container.end());
		m_Size = _Newsize;
		m_MutexLock.Unlock();
	}

	void resize(size_type _Newsize, TElement _Val) {	// determine new length, padding with _Val elements as needed
		m_MutexLock.Lock();
		m_Size = _Newsize;
		m_Container.clear;
		m_Containter[ 0 ] = _Val;
		m_MutexLock.Unlock();
	}

	size_type size() const	{	// return length of sequence
		return m_Size;
	}
/*
	size_type max_size() const
		{	// return maximum possible length of sequence
		return (this->_Alval.max_size());
	}
*/
	bool empty() const	{	// test if sequence is empty
		return (size() == 0);
	}

	const_reference at(size_type _Pos) const {	// subscript nonmutable sequence with checking
		if (size() <= _Pos)
			_Xran();
		return (*(begin() + _Pos));
	}

	reference at(size_type _Pos) {	// subscript mutable sequence with checking
		if (size() <= _Pos)
			_Xran();
		return (*(begin() + _Pos));
	}

	const_reference Get( TIndex ind ) const {
		m_MutexLock.Lock();
		ContainerType::const_iterator it = m_Container.lower_bound( ind );
		if (it->first == ind) {
			m_MutexLock.Unlock();
			return it->second;
		} else {
			if (it != m_Container.begin() ) --it;
			m_MutexLock.Unlock();
			return it->second;
		}
	}

	void Set( const TElement &val, TIndex ind) {
	int range = 3;
	TIndex start = std::max( ind-range, 0 );
	TIndex end = std::max( ind+range, m_Size );
	std::vector< TElement > cont1( begin()+start, begin()+end );
		m_MutexLock.Lock();
		_Set( val, ind );
		m_MutexLock.Unlock();
	TElement ov = cont1[ range ];
	cont1[range] = val;
	std::vector< TElement > cont2( begin()+start, begin()+end );
	if (cont1 != cont2) {
		cont1[range] = ov;
		std::cerr << "setting to:" << val << std::endl;
		for(int i = 0; i < range*2 +1; ++i) {
			std::cerr << cont1[i] << "|" << cont2[i] << "   ";
		}
		std::cerr << std::endl;
		exit(0);
	} 
	}

	void _Set( const TElement &val, TIndex ind) {
		ContainerType::iterator it = m_Container.lower_bound( ind );
		if (it->first == ind) {
			if (it->second == val) return;
			else {
				if (ind < m_Size-1) {
					ContainerType::iterator itPost = it; ++itPost;
					if (itPost->first == ind+1) {
						if (it != m_Container.begin() ) {
							ContainerType::const_iterator itPre = it; --itPre;
							if (itPre->second == val) m_Container.erase( it );
							else it->second = val;
						} else	it->second = val;
						if ( itPost->second == val)
							m_Container.erase( itPost );
					} else {
						TElement tval = it->second;
						if (it != m_Container.begin() ) {
							ContainerType::const_iterator itPre = it; --itPre;
							if (itPre->second == val) m_Container.erase( it );
							else it->second = val;
						} else	it->second = val;
						m_Container[ ind+1 ] = tval;
					}
				} else {
					if (it != m_Container.begin() ) {
						ContainerType::const_iterator itPre = it; --itPre;
						if (itPre->second == val) m_Container.erase( it );
						else it->second = val;
					} else	it->second = val;
				}
			}
		} else {
			ContainerType::const_iterator itPre = it; --itPre;
			if (itPre->second == val) return;
			else {
				if (it->first == ind+1) {
					m_Container[ ind ] = val;
					if ( it->second == val)
						m_Container.erase( it );
				} else {
					m_Container[ ind ] = val;
					m_Container[ ind+1 ] = itPre->second;
				}
			}
		}
	}

	const_reference operator[](size_type _Pos) const {	// subscript nonmutable sequence
		return Get( _Pos );
	}

	reference operator[](size_type _Pos) {	// subscript mutable sequence
		return (*(begin() + _Pos));
	}

	reference front() {	// return first element of mutable sequence
		return (*begin());
	}

	const_reference front() const {	// return first element of nonmutable sequence
		return (*begin());
	}

	reference back() {	// return last element of mutable sequence
		return (*(end() - 1));
	}

	const_reference back() const {	// return last element of nonmutable sequence
		return (*(end() - 1));
	}
/*
	void push_back(const _Ty& _Val)
		{	// insert element at end
		if (size() < capacity())
			_Mylast = _Ufill(_Mylast, 1, _Val);
		else
			insert(end(), _Val);
		}

	void pop_back()
		{	// erase element at end
		if (!empty())
			{	// erase last element
			_Destroy(_Mylast - 1, _Mylast);
			--_Mylast;
			}
		}
*/
/*
	template<class _Iter>
		void assign(_Iter _First, _Iter _Last)
		{	// assign [_First, _Last)
		_Assign(_First, _Last, _Iter_cat(_First));
		}

	template<class _Iter>
		void _Assign(_Iter _Count, _Iter _Val, _Int_iterator_tag)
		{	// assign _Count * _Val
		_Assign_n((size_type)_Count, (_Ty)_Val);
		}

	template<class _Iter>
		void _Assign(_Iter _First, _Iter _Last, input_iterator_tag)
		{	// assign [_First, _Last), input iterators
		erase(begin(), end());
		insert(begin(), _First, _Last);
		}

	void assign(size_type _Count, const _Ty& _Val)
		{	// assign _Count * _Val
		_Assign_n(_Count, _Val);
		}

	iterator insert(iterator _Where, const _Ty& _Val)
		{	// insert _Val at _Where
		size_type _Off = size() == 0 ? 0 : _Where - begin();
		_Insert_n(_Where, (size_type)1, _Val);
		return (begin() + _Off);
		}

	void insert(iterator _Where, size_type _Count, const _Ty& _Val)
		{	// insert _Count * _Val at _Where
		_Insert_n(_Where, _Count, _Val);
		}

	template<class _Iter>
		void insert(iterator _Where, _Iter _First, _Iter _Last)
		{	// insert [_First, _Last) at _Where
		_Insert(_Where, _First, _Last, _Iter_cat(_First));
		}

	template<class _Iter>
		void _Insert(iterator _Where, _Iter _First, _Iter _Last,
			_Int_iterator_tag)
		{	// insert _Count * _Val at _Where
		_Insert_n(_Where, (size_type)_First, (_Ty)_Last);
		}

	template<class _Iter>
		void _Insert(iterator _Where, _Iter _First, _Iter _Last,
			input_iterator_tag)
		{	// insert [_First, _Last) at _Where, input iterators
		for (; _First != _Last; ++_First, ++_Where)
			_Where = insert(_Where, *_First);
		}

	template<class _Iter>
		void _Insert(iterator _Where,
			_Iter _First, _Iter _Last, forward_iterator_tag)
		{	// insert [_First, _Last) at _Where, forward iterators
		size_type _Count = 0;
		_Distance(_First, _Last, _Count);
		size_type _Capacity = capacity();

		if (_Count == 0)
			;
		else if (max_size() - size() < _Count)
			_Xlen();	// result too long
		else if (_Capacity < size() + _Count)
			{	// not enough room, reallocate
			_Capacity = max_size() - _Capacity / 2 < _Capacity
				? 0 : _Capacity + _Capacity / 2;	// try to grow by 50%
			if (_Capacity < size() + _Count)
				_Capacity = size() + _Count;
			pointer _Newvec = this->_Alval.allocate(_Capacity);
			pointer _Ptr = _Newvec;

			_TRY_BEGIN
			_Ptr = _Ucopy(_Myfirst, _ITER_BASE(_Where),
				_Newvec);	// copy prefix
			_Ptr = _Ucopy(_First, _Last, _Ptr);	// add new stuff
			_Ucopy(_ITER_BASE(_Where), _Mylast, _Ptr);	// copy suffix
			_CATCH_ALL
			_Destroy(_Newvec, _Ptr);
			this->_Alval.deallocate(_Newvec, _Capacity);
			_RERAISE;
			_CATCH_END

			_Count += size();
			if (_Myfirst != 0)
				{	// destroy and deallocate old array
				_Destroy(_Myfirst, _Mylast);
				this->_Alval.deallocate(_Myfirst, _Myend - _Myfirst);
				}
			_Myend = _Newvec + _Capacity;
			_Mylast = _Newvec + _Count;
			_Myfirst = _Newvec;
			}
		else if ((size_type)(end() - _Where) < _Count)
			{	// new stuff spills off end
			_Ucopy(_ITER_BASE(_Where), _Mylast,
				_ITER_BASE(_Where) + _Count);	// copy suffix
			_Iter _Mid = _First;
			advance(_Mid, end() - _Where);

			_TRY_BEGIN
			_Ucopy(_Mid, _Last, _Mylast);	// insert new stuff off end
			_CATCH_ALL
			_Destroy(_ITER_BASE(_Where) + _Count, _Mylast + _Count);
			_RERAISE;
			_CATCH_END

			_Mylast += _Count;
			copy(_First, _Mid, _ITER_BASE(_Where));	// insert to old end
			}
		else
			{	// new stuff can all be assigned
			pointer _Oldend = _Mylast;
			_Mylast = _Ucopy(_Oldend - _Count, _Oldend,
				_Mylast);	// copy suffix
			copy_backward(_ITER_BASE(_Where), _Oldend - _Count,
				_Oldend);	// copy hole
			copy(_First, _Last, _ITER_BASE(_Where));	// insert into hole
			}
		}

	iterator erase(iterator _Where)
		{	// erase element at where
		copy(_ITER_BASE(_Where) + 1, _Mylast, _ITER_BASE(_Where));
		_Destroy(_Mylast - 1, _Mylast);
		--_Mylast;
		return (_Where);
		}

	iterator erase(iterator _First, iterator _Last)
		{	// erase [_First, _Last)
		if (_First != _Last)
			{	// worth doing, copy down over hole
			pointer _Ptr = copy(_ITER_BASE(_Last), _Mylast,
				_ITER_BASE(_First));
			_Destroy(_Ptr, _Mylast);
			_Mylast = _Ptr;
			}
		return (_First);
		}
*/
	void clear() {	// erase all
		m_MutexLock.Lock();
		m_Size = 0;
		m_Container.clear();
		m_MutexLock.Unlock();
	}

	void swap(sparsefield& _Right)	{	// exchange contents with _Right
		m_MutexLock.Lock();
		std::swap(m_Size, _Right.m_Size);
		std::swap(m_Container, _Right.m_Container);
		m_MutexLock.Unlock();
	}
	

	void Print(std::ostream & out) {
		out << "sparsefield::size = " << m_Size << std::endl;
		for(sparsefield::ContainerType::const_iterator ci = m_Container.begin(); ci!=m_Container.end(); ++ci) {
			out << "[" << ci->first << "]->" << ci->second << "  ";
		}
		out << std::endl;
	}
	
protected:
	TIndex m_Size;
	ContainerType m_Container;
	itk::SimpleFastMutexLock m_MutexLock;
};


	


#endif 

