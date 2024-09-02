#include <cstdint>
#include <cassert>
#include <cstdio>
#include <array>
#include <algorithm>
#include <limits>


template<size_t N>
class TDigest {
	static_assert(N >= 2, "if N is less 2, there will be nothing to compress to.");

	struct Centroid {
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
			if (*this)
				printf("Addr %p mean: %10.4f | weight: %5zu\n", (void *) this, mean, weight);
			else
				printf("---empty---\n");
		}

		friend bool operator<(Centroid const &a, Centroid const &b){
			return a.mean < b.mean;
		}
	};

	using CentroidArray = std::array<Centroid, N>;

	CentroidArray	centroids_;
	double		delta_;

public:
	constexpr TDigest(double delta) : delta_(delta){}

	constexpr static size_t size(){
		return N;
	}

	enum class Compression{
		NONE		,
		STANDARD	,
		AGGRESSIVE
	};

	template<Compression C = Compression::AGGRESSIVE>
	void add(double value, uint64_t weight = 1){
		assert(weight > 0);

		constexpr bool update_follows = true;

		auto update = [&](size_t count){
			if (count < size())
				centroids_[count] = { value, weight };

			if (count + 1 < size())
				centroids_[count + 1].clear();
		};

		auto count = getCount_();

		if (count < size())
			return update(count);

		if constexpr(C == Compression::NONE)
			return;

		if constexpr(C == Compression::STANDARD)
			count = compressNormal_<update_follows>(count);

		if constexpr(C == Compression::AGGRESSIVE)
			count = compressAggressive_<update_follows>(count);

		// drop the value
		// should never happen, if Aggressive,

		return update(count);
	}

	void compress(){
		constexpr bool update_follows = true;

		auto const count = getCount_();

		[[maybe_unused]]
		auto const x = compressNormal_<not update_follows>(count);
	}

	void print() const{
		printf("Centroids:\n");

		for (auto &x : centroids_)
			if (x)
				x.print();
			else
				break;
	}

	void printAll() const{
		printf("Centroids:\n");

		for (auto &x : centroids_)
			x.print();
	}


private:
	[[nodiscard]]
	size_t getCount_() const{
		for(size_t i = 0; i < size(); ++i)
			if (!centroids_[i])
				return i;

		return size();
	}

	template<bool AddFollows>
	[[nodiscard]]
	size_t compressNormal_(size_t count){
		if (count < 2)
			return count;

		sort_(count);

		return compressCentroids_<1,AddFollows>(count, delta_);
	}

	template<bool AddFollows>
	[[nodiscard]]
	size_t compressAggressive_(size_t count){
		if (count < 2)
			return count;

		sort_(count);

		auto const distance = findMinDistance_(count);

		if (distance > delta_)
			return compressCentroids_<0,AddFollows>(count, distance);
		else
			return compressCentroids_<1,AddFollows>(count, delta_);
	}

	void sort_(size_t count){
		std::sort(std::begin(centroids_), std::begin(centroids_) + count);
	}

	[[nodiscard]]
	double findMinDistance_(size_t count) const{
		assert(count > 1);

		double minDistance = std::numeric_limits<double>::max();

		for (size_t i = 0; i < count - 1; ++i){
			auto const distance = std::abs(centroids_[i].mean - centroids_[i + 1].mean);

			if (distance < minDistance)
				minDistance = distance;
		}

		return minDistance;
	}

	template<bool UseWeight, bool AddFollows>
	[[nodiscard]]
	size_t compressCentroids_(size_t count, double delta){
		assert(count > 1);

		if constexpr(!UseWeight){
		//	printf("compressAggressive with delta: %10.5f\n", delta);
		}

	//	printf("delta: %10.5f\n", delta);

		size_t newCount = 0;
		auto   current = centroids_[0];

		auto const _ = [](double weight) -> double{
			if constexpr(UseWeight)
				return weight;
			else
				return 1.0;
		};

		for (size_t i = 1; i < count; ++i){
			auto const distance = std::abs(centroids_[i].mean - current.mean);
			auto const weight_u = current.weight + centroids_[i].weight;
			auto const weight   = static_cast<double>(weight_u);

			if (_(weight) * distance <= delta) {
				current.mean   = (current.getWeightedMean() + centroids_[i].getWeightedMean()) / weight;
				current.weight = weight_u;
			}else{
				centroids_[newCount++] = current;
				current = centroids_[i];
			}
		}

		centroids_[newCount++] = current;

		if (!AddFollows){
			// now we need to invalidate *ONE* centroid after newCount
			// because there *WILL* be insert coming

			if (newCount < size())
				centroids_[newCount].clear();
		}else{
			// insert will invalidate as needed
		}

		return newCount;
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
	printf("View...\n");
	auto td = get();
	td.printAll();



	printf("None...\n");
	td = get();
	td.add<C::NONE>(1.52);
	td.printAll();



	printf("Standard...\n");
	td = get();
	td.add<C::STANDARD>(1.52);
	td.printAll();



	printf("Aggresive...\n");
	td = get();
	td.add<C::AGGRESSIVE>(1.52);
	td.printAll();



	printf("Bad input + standard...\n");
	td = getBad();
	td.add<C::STANDARD>(1.52);
	td.printAll();



	printf("Bad input + aggresive...\n");
	td = getBad();
	td.add<C::AGGRESSIVE>(1.52);
	td.printAll();



	printf("Extreme input...\n");
	td = getExtreme();
	td.compress();
	td.printAll();



	printf("Try to smash the array...\n");
	td = getExtreme();
	td.compress();
	td.add(1.5);
	td.printAll();
}

