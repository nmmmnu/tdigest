#include "tdigest.h"

#include <cstdio>
#include <iterator>

namespace{
	constexpr size_t SIZE  = 5;
	constexpr double DELTA = 0.05;

	constexpr auto CM = RawTDigest::Compression::NONE;

	using Centroid = RawTDigest::Centroid;

	auto get(Centroid *cd){
		RawTDigest td{ SIZE, DELTA };

		td.clearFast(cd);

		td.add<CM>(cd, 25.0);
		td.add<CM>(cd, 1.50);
		td.add<CM>(cd, 1.51);
		td.add<CM>(cd, 1.46);
		td.add<CM>(cd, 1.47);
	}

	auto getBad(Centroid *cd){
		RawTDigest td{ SIZE, DELTA };

		td.clearFast(cd);

		td.add<CM>(cd, 31);
		td.add<CM>(cd, 41);
		td.add<CM>(cd, 51);
		td.add<CM>(cd, 10);
		td.add<CM>(cd, 20);
	}

	auto getExtreme(Centroid *cd){
		RawTDigest td{ SIZE, DELTA };

		td.clearFast(cd);

		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
	}

} // anonymous namespace

int main(){
	RawTDigest td{ SIZE, DELTA };

	Centroid *cd = reinterpret_cast<Centroid *>(malloc(td.bytes()));

	using C = RawTDigest::Compression;

	constexpr double pp[3]{ 0.05, 0.50, 0.95 };
	double oo[3];



	printf("None...\n");
	get(cd);
	td.add<C::NONE>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Standard...\n");
	get(cd);
	td.add<C::STANDARD>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Aggresive...\n");
	get(cd);
	td.add<C::AGGRESSIVE>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Bad input + standard...\n");
	getBad(cd);
	td.add<C::STANDARD>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));


	printf("Bad input + aggresive...\n");
	getBad(cd);
	td.add<C::AGGRESSIVE>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));

	td.percentile(cd, std::cbegin(pp), std::cend(pp), std::begin(oo));
	for(auto const &x : oo)
		printf("-> %10.6f\n", x);



	printf("Extreme input...\n");
	getExtreme(cd);
	td.compress(cd);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));

	td.percentile(cd, std::cbegin(pp), std::cend(pp), std::begin(oo));
	for(auto const &x : oo)
		printf("-> %10.6f\n", x);



	free(cd);
}


