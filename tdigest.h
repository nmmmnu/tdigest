#ifndef T_DIGEST_H_
#define T_DIGEST_H_

#include <cstdint>
#include <cassert>
#include <cstring>
#include <algorithm>	// transform

class RawTDigest{
	size_t	capacity_;
	double	delta_;

	static const size_t sizeof_Centroid__;

public:
	struct Centroid;

public:
	constexpr RawTDigest(size_t capacity, double delta) : capacity_(capacity), delta_(delta){
		assert(capacity_ >= 2);
	}

	enum class Compression{
		NONE		,
		STANDARD	,
		AGGRESSIVE
	};

	constexpr size_t capacity() const{
		return capacity_;
	}

	constexpr size_t bytes() const{
		return capacity_ * sizeof_Centroid__;
	}

	void print(const Centroid *cd) const;

public:
	static void clearFast(Centroid *cd){
		memset(cd, 0, sizeof_Centroid__);
	}

	void clear(Centroid *cd) const{
		memset(cd, 0, bytes());
	}

	void load(Centroid *cd, const void *src) const{
		memcpy(cd, src, bytes());
	}

	void store(const Centroid *cd, void *dest) const{
		memcpy(dest, cd, bytes());
	}

public:
	template<Compression C = Compression::AGGRESSIVE>
	void add(Centroid *cd, double value, uint64_t weight = 1) const;

	size_t compress(Centroid *cd) const{
		size_t const size = getSize_(cd);

		return compressNormal_(cd, size);
	}

	double percentile_50(const Centroid *cd) const{
		return percentile(cd, 0.50);
	}

	double percentile_95(const Centroid *cd) const{
		return percentile(cd, 0.95);
	}

	double percentile(const Centroid *cd, double const p) const{
		assert(p >= 0.00 && p <= 1.00);

		auto [weight, size] = getWeightAndSize_(cd);

		return percentile_(cd, size, weight, p);
	}

	template<typename IT, typename OutIT>
	void percentile(const Centroid *cd, IT first, IT last, OutIT out) const{
		auto [weight, size] = getWeightAndSize_(cd);

		auto f = [&](double p){
			assert(p >= 0.00 && p <= 1.00);
			return percentile_(cd, size, weight, p);
		};

		std::transform(first, last, out, f);
	}

private:
	size_t getSize_(const Centroid *cd) const;

	std::pair<uint64_t, size_t> getWeightAndSize_(const Centroid *cd) const;

	double percentile_(const Centroid *cd, size_t size, uint64_t weight, double const p) const;

	size_t compressNormal_(Centroid *cd, size_t size) const;

	size_t compressAggressive_(Centroid *cd, size_t size) const;

	template<bool UseWeight>
	size_t compressCentroids_(Centroid *cd, size_t size, double delta) const;

	static double findMinDistance__(const Centroid *cd, size_t size);
};

#endif

