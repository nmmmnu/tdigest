#include <cstdint>
#include <cassert>
#include <cstdio>

#include <algorithm>
#include <iterator>
#include <limits>



namespace t_digest_impl_{
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



template<size_t N>
class TDigest {
	static_assert(N >= 2, "if N is less 2, there will be nothing to compress to.");

	struct Centroid{
		double   mean	= 0;
		uint64_t weight	= 0;

		constexpr Centroid() = default;
		constexpr Centroid(double m, uint64_t w) : mean(m), weight(w){}

		void clear(){
			mean   = 0;
			weight = 0;
		}

		operator bool() const{
			return weight;
		}

		double getWeightedMean() const{
			return mean * static_cast<double>(weight);
		}

		void print() const{
			printf("> Addr %p | mean: %10.4f | weight: %5zu\n", (void *) this, mean, weight);
		}

		friend bool operator<(Centroid const &a, Centroid const &b){
			return a.mean < b.mean;
		}
	};

	double		delta_;

	uint64_t	size_		= 0;
	uint64_t	weight_		= 0;

	Centroid	centroids_[N];

private:
	constexpr auto begin(){
		return std::begin(centroids_);
	}

	constexpr auto end(){
		return std::begin(centroids_) + size();
	}

	constexpr auto begin() const{
		return std::begin(centroids_);
	}

	constexpr auto end() const{
		return std::begin(centroids_) + size();
	}

public:
	constexpr TDigest(double delta) : delta_(delta){}

	constexpr static size_t capacity(){
		return N;
	}

	constexpr uint64_t size() const{
		return size_;
	}

	enum class Compression{
		NONE		,
		STANDARD	,
		AGGRESSIVE
	};

	template<Compression C = Compression::AGGRESSIVE>
	void add(double value, uint64_t weight = 1){
		assert(weight > 0);

		auto update = [&](){
			using namespace t_digest_impl_;

			insertIntoSortedRange(begin(), end(), Centroid{ value, weight });

			++size_;

			weight_ += weight;
		};

		if (size() < capacity())
			return update();

		if constexpr(C == Compression::NONE)
			return;

		if constexpr(C == Compression::STANDARD)
			compressNormal_();

		if constexpr(C == Compression::AGGRESSIVE)
			compressAggressive_();

		if (size() < capacity())
			return update();

		// drop the value
		// should be unreachible if Aggressive,
	}

	void compress(){
		return compressNormal_();
	}

	void print() const{
		printf("Centroids, size %zu, capacity %zu\n", size(), capacity());

		for(auto it = begin(); it != end(); ++it)
			it->print();
	}

	double percentile_50() const{
		return percentile(0.50);
	}

	double percentile_95() const{
		return percentile(0.95);
	}

	double percentile(double const p) const{
		assert(p >= 0.00 && p <= 1.00);

		if (size() == 0)
			return 0;

		double const targetRank = p * static_cast<double>(weight_);
		double cumulativeWeight = 0;

		for(auto it = begin(); it != end(); ++it){
			cumulativeWeight += static_cast<double>(it->weight);

			if (cumulativeWeight > targetRank)
				return it->mean;
		}

		return std::prev(end())->mean;
	}

private:
	void compressNormal_(){
		if (size() < 2)
			return;

		return compressCentroids_<1>(delta_);
	}

	void compressAggressive_(){
		if (size() < 2)
			return;

		auto const distance = findMinDistance_();

		if (distance > delta_)
			return compressCentroids_<0>(distance);
		else
			return compressCentroids_<1>(delta_);
	}

	double findMinDistance_() const{
		double minDistance = std::numeric_limits<double>::max();

		for(auto it = begin(); it != end(); ++it){
			auto const distance = std::abs(it->mean - std::next(it)->mean);

			if (distance < minDistance)
				minDistance = distance;
		}

		return minDistance;
	}

	template<bool UseWeight>
	void compressCentroids_(double delta){
		assert(size() > 1);

		if constexpr(!UseWeight){
		//	printf("compressAggressive with delta: %10.5f\n", delta);
		}

		size_t newSize = 0;
		auto   current = centroids_[0];

		auto const _ = [](double weight) -> double{
			if constexpr(UseWeight)
				return weight;
			else
				return 1.0;
		};

		for (size_t i = 1; i < size(); ++i){
			auto const distance = std::abs(centroids_[i].mean - current.mean);
			auto const weight_u = current.weight + centroids_[i].weight;
			auto const weight   = static_cast<double>(weight_u);

			if (_(weight) * distance <= delta) {
				current.mean   = (current.getWeightedMean() + centroids_[i].getWeightedMean()) / weight;
				current.weight = weight_u;
			}else{
				centroids_[newSize++] = current;
				current = centroids_[i];
			}
		}

		centroids_[newSize++] = current;

		size_		= newSize;
	}
};



namespace{
	using MyTDigest = TDigest<5>;
	using C = MyTDigest::Compression;

	constexpr double DELTA = 0.05;

	auto get(double delta = DELTA){
		MyTDigest td{ delta };

		td.add<C::NONE>(25.0);

		td.add<C::NONE>(1.50);
		td.add<C::NONE>(1.51);

		td.add<C::NONE>(1.46);
		td.add<C::NONE>(1.47);

		return td;
	}

	auto getBad(double delta = DELTA){
		MyTDigest td{ delta };

		td.add<C::NONE>(31);
		td.add<C::NONE>(41);
		td.add<C::NONE>(51);

		td.add<C::NONE>(10);
		td.add<C::NONE>(20);

		return td;
	}

	auto getExtreme(double delta = DELTA){
		MyTDigest td{ delta };

		td.add<C::NONE>(1.5);
		td.add<C::NONE>(1.5);
		td.add<C::NONE>(1.5);
		td.add<C::NONE>(1.5);
		td.add<C::NONE>(1.5);

		return td;
	}

} // anonymous namespace

int main() {
	printf("None...\n");
	auto td = get();
	td.add<C::NONE>(1.52);
	td.print();
	printf("%10.6f\n", td.percentile_50());



	printf("Standard...\n");
	td = get();
	td.add<C::STANDARD>(1.52);
	td.print();
	printf("%10.6f\n", td.percentile_50());



	printf("Aggresive...\n");
	td = get();
	td.add<C::AGGRESSIVE>(1.52);
	td.print();
	printf("%10.6f\n", td.percentile_50());



	printf("Bad input + standard...\n");
	td = getBad();
	td.add<C::STANDARD>(1.52);
	td.print();
	printf("%10.6f\n", td.percentile_50());



	printf("Bad input + aggresive...\n");
	td = getBad();
	td.add<C::AGGRESSIVE>(1.52);
	td.print();
	printf("%10.6f\n", td.percentile_50());



	printf("Extreme input...\n");
	td = getExtreme();
	td.compress();
	td.print();
	printf("%10.6f\n", td.percentile_50());
}


