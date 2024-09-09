#include "tdigest.h"

#include <limits>
#include <cstdio>

namespace {
	#if 0
	template<typename T>
	void insertSorted2(T *arr, size_t const n, T &&last){
	    auto it = std::lower_bound(arr, arr + n, last);
	    std::move_backward(it, arr + n, arr + n + 1);
	    *it = std::move(last);
	}
	#endif

	template<typename IT>
	void insertIntoSortedRange(IT first, IT last, typename std::iterator_traits<IT>::value_type &&item){
	    auto it = std::lower_bound(first, last, item);
	    std::move_backward(it, last, std::next(last));

	    *it = std::move(item);
	}
}



struct RawTDigest::Centroid{
	double   mean_;
	uint64_t weight_;

	constexpr static auto create(double mean, uint64_t weight){
		return Centroid{ mean, weight };
	}

	constexpr void clear(){
		mean_   = 0;
		weight_ = 0;
	}

	constexpr auto getMean() const{
		return mean_;
	}

	constexpr auto getWeight() const{
		return weight_;
	}

	constexpr operator bool() const{
		return weight_;
	}

	constexpr double getWeightedMean() const{
		return getMean() * static_cast<double>(getWeight());
	}

	void print() const{
		printf("> Addr %p | mean: %10.4f | weight: %5zu\n", (void *) this, getMean(), getWeight());
	}

	friend constexpr bool operator<(Centroid const &a, Centroid const &b){
		return a.getMean() < b.getMean();
	}
};

static_assert(std::is_trivial_v<RawTDigest::Centroid>);

const size_t RawTDigest::sizeof_Centroid__ = sizeof(RawTDigest::Centroid);



void RawTDigest::print(const Centroid *cd) const{
	printf("Centroids, capacity %zu\n", capacity());

	for(size_t i = 0; i < capacity(); ++i){
		auto const &x = cd[i];
		if (!x)
			break;

		x.print();
	}
}



size_t RawTDigest::getSize_(const Centroid *cd) const{
	size_t size = 0;

	for(size_t i = 0; i < capacity(); ++i){
		auto const &x = cd[i];
		if (!x)
			break;

		++size;
	}

	return size;
}

std::pair<uint64_t, size_t> RawTDigest::getWeightAndSize_(const Centroid *cd) const{
	std::pair<uint64_t, size_t> r{ 0, 0 };

	for(size_t i = 0; i < capacity(); ++i){
		auto const &x = cd[i];
		if (!x)
			break;

		r.first += x.getWeight();
		++r.second;
	}

	return r;
}

double RawTDigest::percentile_(const Centroid *cd, size_t size, uint64_t weight, double const p) const{
	if (size == 0)
		return 0;

	double const targetRank = p * static_cast<double>(weight);
	double cumulativeWeight = 0;

	for(size_t i = 0; i < size - 1; ++i){
		cumulativeWeight += static_cast<double>(cd[i].getWeight());
		if (cumulativeWeight >= targetRank)
			return cd[i].getMean();
	}

	return cd[size - 1].getMean();
}

template<RawTDigest::Compression C>
void RawTDigest::add(Centroid *cd, double value, uint64_t weight) const{
	assert(weight > 0);

	auto size = getSize_(cd);

	auto insert = [&](){
		insertIntoSortedRange(cd, cd + size, Centroid::create(value, weight) );

		if (++size < capacity())
			cd[size].clear();
	};

	if (size < capacity())
		return insert();

	if constexpr(C == Compression::NONE)
		return;

	if constexpr(C == Compression::STANDARD)
		size = compressNormal_(cd, size);

	if constexpr(C == Compression::AGGRESSIVE)
		size = compressAggressive_(cd, size);

	if (size < capacity())
		return insert();

	// drop the value
	// should be unreachible if Aggressive,
}

template void RawTDigest::add<RawTDigest::Compression::NONE		>(Centroid *cd, double value, uint64_t weight) const;
template void RawTDigest::add<RawTDigest::Compression::STANDARD		>(Centroid *cd, double value, uint64_t weight) const;
template void RawTDigest::add<RawTDigest::Compression::AGGRESSIVE	>(Centroid *cd, double value, uint64_t weight) const;



size_t RawTDigest::compressNormal_(Centroid *cd, size_t size) const{
	if (size < 2)
		return size;

	return compressCentroids_<1>(cd, size, delta_);
}

size_t RawTDigest::compressAggressive_(Centroid *cd, size_t size) const{
	if (size < 2)
		return size;

	auto const distance = findMinDistance__(cd, size);

	if (distance > delta_)
		return compressCentroids_<0>(cd, size, distance);
	else
		return compressCentroids_<1>(cd, size, delta_);
}

template<bool UseWeight>
size_t RawTDigest::compressCentroids_(Centroid *cd, size_t size, double delta) const{
	assert(size > 1);

	size_t newSize = 0;
	auto   current = cd[0];

	auto const _ = [](double weight) -> double{
		if constexpr(UseWeight)
			return weight;
		else
			return 1.0;
	};

	for (size_t i = 1; i < size; ++i){
		auto const distance = std::abs(cd[i].getMean() - current.getMean());
		auto const weight_u = current.getWeight() + cd[i].getWeight();
		auto const weight   = static_cast<double>(weight_u);

		if (_(weight) * distance <= delta) {
			current = Centroid::create(
					(current.getWeightedMean() + cd[i].getWeightedMean()) / weight,
					weight_u
			);
		}else{
			cd[newSize++] = current;
			current = cd[i];
		}
	}

	cd[newSize++] = current;

	if (newSize < capacity())
		cd[newSize].clear();

	return newSize;
}



double RawTDigest::findMinDistance__(const Centroid *cd, size_t const size){
	assert(size > 1);

	double minDistance = std::numeric_limits<double>::max();

	for(auto it = cd; it != cd + size - 1; ++it){
		auto const distance = std::abs(it->getMean() - std::next(it)->getMean());

		if (distance < minDistance)
			minDistance = distance;
	}

	return minDistance;
}



#if 0

class MultiPercentile{
	RawTDigest const	&td;
	const Centroid		*cd;
	uint64_t		weight;
	size_t			size;

public:
	MultiPercentile(RawTDigest const &td, const Centroid *cd);
	double operator()(double p) const;
};


RawTDigest::MultiPercentile::MultiPercentile(RawTDigest const &td, const Centroid *cd) : td(td), cd(cd){
	auto [weight, size] = td.getWeightAndSize_(cd);
	this->weight = weight;
	this->size   = size;
}

double RawTDigest::MultiPercentile::operator()(double p) const{
	assert(p >= 0.00 && p <= 1.00);
	return td.percentile_(cd, size, weight, p);
}

#endif

