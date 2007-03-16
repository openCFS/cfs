#ifndef _POLYGONITERATOR_HH_2007_
#define _POLYGONITERATOR_HH_2007_

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"

namespace CoupledField
{
    // define a type for polygons
    typedef StdVector< Vector<Double> > Polygon;

    class PolygonIterator
    {

	public:

	    inline PolygonIterator(Polygon &p)
		: poly_(&p), pos_(0), start_(0), dir_(1), dirty_(false) {};

	    inline PolygonIterator(Polygon &p, UInt pos)
		: poly_(&p), dir_(1), dirty_(false)
	    {
		if (pos < p.GetSize())
		    pos_ = pos;
		else
		    pos_ = 0;
		start_ = pos_;
	    };
	    
	    inline ~PolygonIterator() {};

	    inline void Swap(PolygonIterator &pi)
	    {
		Polygon *tmp_poly = poly_;
		poly_ = pi.poly_;
		pi.poly_ = tmp_poly;

		UInt tmp = pos_;
		pos_ = pi.pos_;
		pi.pos_ = tmp;

		tmp = start_;
		start_ = pi.start_;
		pi.start_ = tmp;

		tmp = dir_;
		dir_ = pi.dir_;
		pi.dir_ = tmp;

		tmp = (UInt) dirty_;
		dirty_ = pi.dirty_;
		pi.dirty_ = (bool) tmp;
	    };
	    
	    inline void Reverse() { dir_ *= -1; };

	    inline Polygon& GetPolygon() const { return *poly_; };

	    inline UInt GetPos() const { return pos_; };
	    
	    inline bool AtBegin() const { return (pos_ == start_); };
	    
	    inline bool AtEnd() const
	    {
		return ( (pos_ == start_) && dirty_ );
	    };
	    
	    inline void SetBegin(UInt b)
	    {
		if ((b >= 0) && (b < poly_->GetSize()))
		    start_ = b;
	    };
	    
	    inline void SetBegin() { start_ = pos_; };
	    
	    inline void Seek(UInt pos)
	    {
		if ((pos >= 0) && (pos < poly_->GetSize()))
		    pos_ = pos;
		dirty_ = false;
	    };
	    
	    inline void Rewind() { pos_ = start_; dirty_ = false; };
	    
	    inline Vector<Double>& Prev(UInt dec = 1) const
	    {
		Integer pos = pos_ - dir_ * dec;
		if (pos >= (Integer) poly_->GetSize())
		    pos %= poly_->GetSize();
		else
		    while (pos < 0)
			pos += poly_->GetSize();
		return (*poly_)[pos];
	    };
	    
	    inline Vector<Double>& Next(UInt inc = 1) const
	    {
		Integer pos = pos_ + dir_ * inc;
		if (pos >= (Integer) poly_->GetSize())
		    pos %= poly_->GetSize();
		else
		    while (pos < 0)
			pos += poly_->GetSize();
		return (*poly_)[pos];
	    };

	    inline PolygonIterator& operator--()
	    {
		Integer pos = pos_ - dir_;
		if (pos < 0)
		    pos += poly_->GetSize();
		else if (pos >= (Integer) poly_->GetSize())
		    pos -= poly_->GetSize();
		pos_ = (UInt) pos;
		dirty_ = true;
		return *this;
	    };
	    
	    inline PolygonIterator& operator++()
	    {
		Integer pos = pos_ + dir_;
		if (pos >= (Integer) poly_->GetSize())
		    pos -= poly_->GetSize();
		else if (pos < 0)
		    pos += poly_->GetSize();
		pos_ = (UInt) pos;
		dirty_ = true;
		return *this;
	    };
	    
	    inline Vector<Double>& operator*() const
	    {
		return (*poly_)[pos_];
	    };

	private:
	    
	    Polygon *poly_;
	    UInt pos_;
	    UInt start_;
	    Integer dir_;
	    bool dirty_;
	    
    }; // class PolygonIterator

    class ConstPolygonIterator {
	
	public:
	    
	    inline ConstPolygonIterator(const Polygon &p)
		: pos_(0), start_(0), dir_(1), dirty_(false)
	    {
		poly_ = const_cast<Polygon*>(&p);
	    };
	    
	    inline ConstPolygonIterator(const Polygon &p, UInt pos)
		: dir_(1), dirty_(false)
	    {
		if (pos < p.GetSize())
		    pos_ = pos;
		else
		    pos_ = 0;
		start_ = pos_;
		poly_ = const_cast<Polygon*>(&p);
	    };
	    
	    inline ~ConstPolygonIterator() {};

	    inline void Swap(ConstPolygonIterator &pi)
	    {
		Polygon *tmp_poly = poly_;
		poly_ = pi.poly_;
		pi.poly_ = tmp_poly;

		UInt tmp = pos_;
		pos_ = pi.pos_;
		pi.pos_ = tmp;

		tmp = start_;
		start_ = pi.start_;
		pi.start_ = tmp;

		tmp = dir_;
		dir_ = pi.dir_;
		pi.dir_ = tmp;

		tmp = (UInt) dirty_;
		dirty_ = pi.dirty_;
		pi.dirty_ = (bool) tmp;
	    };
	    
	    inline void Reverse() { dir_ *= -1; };

	    inline const Polygon& GetPolygon() const { return *poly_; };

	    inline UInt GetPos() const { return pos_; };
	    
	    inline bool AtBegin() const { return (pos_ == start_); };
	    
	    inline bool AtEnd() const
	    {
		return ( (pos_ == start_) && dirty_ );
	    };
	    
	    inline void SetBegin(UInt b)
	    {
		if ((b >= 0) && (b < poly_->GetSize()))
		    start_ = b;
	    };
	    
	    inline void SetBegin() { start_ = pos_; };
	    
	    inline void Seek(UInt pos)
	    {
		if ((pos >= 0) && (pos < poly_->GetSize()))
		    pos_ = pos;
		dirty_ = false;
	    };
	    
	    inline void Rewind() { pos_ = start_; dirty_ = false; };
	    
	    inline const Vector<Double>& Prev(UInt dec = 1) const
	    {
		Integer pos = pos_ - dir_ * dec;
		if (pos >= (Integer) poly_->GetSize())
		    pos %= poly_->GetSize();
		else
		    while (pos < 0)
			pos += poly_->GetSize();
		return (*poly_)[pos];
	    };
	    
	    inline const Vector<Double>& Next(UInt inc = 1) const
	    {
		Integer pos = pos_ + dir_ * inc;
		if (pos >= (Integer) poly_->GetSize())
		    pos %= poly_->GetSize();
		else
		    while (pos < 0)
			pos += poly_->GetSize();
		return (*poly_)[pos];
	    };

	    inline ConstPolygonIterator& operator--()
	    {
		Integer pos = pos_ - dir_;
		if (pos < 0)
		    pos += poly_->GetSize();
		else if (pos >= (Integer) poly_->GetSize())
		    pos -= poly_->GetSize();
		pos_ = (UInt) pos;
		dirty_ = true;
		return *this;
	    };
	    
	    inline ConstPolygonIterator& operator++()
	    {
		Integer pos = pos_ + dir_;
		if (pos >= (Integer) poly_->GetSize())
		    pos -= poly_->GetSize();
		else if (pos < 0)
		    pos += poly_->GetSize();
		pos_ = (UInt) pos;
		dirty_ = true;
		return *this;
	    };
	    
	    inline const Vector<Double>& operator*() const
	    {
		return (*poly_)[pos_];
	    };

	private:
	    
	    Polygon *poly_;
	    UInt pos_;
	    UInt start_;
	    Integer dir_;
	    bool dirty_;
	    
    }; // class ConstPolygonIterator

} // namespace CoupledField

#endif // _POLYGONITERATOR_HH_2007_
