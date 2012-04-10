// $Id: CoordinateListTest.cpp 3155 2010-12-03 13:57:00Z strk $
// 
// Test Suite for geos::geom::CoordinateList class.

// tut
#include <tut.hpp>
// geos
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateList.h>
#include <geos/geom/CoordinateArraySequence.h>
// std
#include <memory>
#include <string>
#include <vector>

namespace tut
{
    //
    // Test Group
    //

    // Common data used by tests
    struct test_coordinatelist_data
    {

        test_coordinatelist_data() {}
    };

    typedef test_group<test_coordinatelist_data> group;
    typedef group::object object;

    group test_coordinatelist_group("geos::geom::CoordinateList");

    //
    // Test Cases
    //

    // Test insert and erase
    template<>
    template<>
    void object::test<1>()
    {
		using geos::geom::Coordinate;
		
		const Coordinate a(0, 0);
		const Coordinate b(10, 10);
		const Coordinate c(20, 20);
		const Coordinate d(5, 5);

		geos::geom::CoordinateList::iterator it, it2;

		std::auto_ptr< std::vector<Coordinate> > col( new std::vector<Coordinate>() );
		col->push_back(a);
		col->push_back(b);
		col->push_back(c);

		// coordinates are copied
		geos::geom::CoordinateList clist(*col);

		ensure_equals( clist.size(), 3u );

		it = clist.begin();
		clist.insert(++it, d);

		ensure_equals( clist.size(), 4u );

		it = clist.begin();
		ensure_equals(*it, a);
		++it;
		ensure_equals(*it, d);
		++it;
		ensure_equals(*it, b);
		++it;
		ensure_equals(*it, c);


		it = clist.begin();
		++it; ++it;
		clist.erase(it);
		ensure_equals( clist.size(), 3u );

		it = clist.begin();
		ensure_equals(*it, a);
		++it;
		ensure_equals(*it, d);
		++it;
		ensure_equals(*it, c);

		clist.insert(clist.end(), b);

		ensure_equals( clist.size(), 4u );

		it = clist.begin(); ++it;
		it2 = it; ++it2; ++it2;
		clist.erase(it, it2);

		ensure_equals( clist.size(), 2u );

		it = clist.begin();
		ensure_equals(*it, a);
		++it;
		ensure_equals(*it, b);

	}

  // Test insert with and without duplicates
  template<>
  template<>
  void object::test<2>()
  {
    using geos::geom::Coordinate;

    geos::geom::CoordinateList clist;
    ensure_equals( clist.size(), 0u );

    clist.insert(clist.end(), Coordinate(0, 0));
    ensure_equals( clist.size(), 1u );

    clist.insert(clist.end(), Coordinate(0, 0), false);
    ensure_equals( clist.size(), 1u );

    clist.insert(clist.end(), Coordinate(0, 0), true);
    ensure_equals( clist.size(), 2u );

    clist.insert(clist.end(), Coordinate(1, 1), true);
    ensure_equals( clist.size(), 3u );

    geos::geom::CoordinateList::iterator it = clist.end(); 
    --it;
    clist.insert(it, Coordinate(0, 0), false);
    ensure_equals( clist.size(), 3u );
  }


} // namespace tut
