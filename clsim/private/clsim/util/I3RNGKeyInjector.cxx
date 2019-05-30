
#include "icetray/I3ConditionalModule.h"
#include "icetray/I3PODHolder.h"
#include <boost/make_shared.hpp>

class I3RNGKeyInjector : public I3ConditionalModule {
public:
    I3RNGKeyInjector(const I3Context &ctx) : I3ConditionalModule(ctx)
    {
        AddParameter("Dataset", "Dataset number", uint64_t(0));
        AddParameter("Job", "Job number", uint64_t(0));
        AddParameter("Frame", "Starting frame number", uint64_t(0));
        AddParameter("Key", "Name of frame object to create", "RNGKey");
    }
    void Configure()
    {
        uint64_t dataset, job;
        GetParameter("Dataset", dataset);
        GetParameter("Job", job);
        GetParameter("Frame", frameIndex_);
        GetParameter("Key", name_);
        key_ = (dataset << 44) + (job << 24);
    }
    void DAQ(I3FramePtr frame)
    {
        frame->Put(name_, boost::make_shared<I3PODHolder<uint64_t>>(key_ + (frameIndex_ & ((1<<24)-1))));
        frameIndex_++;
        PushFrame(frame);
    }
private:
    std::string name_;
    uint64_t key_;
    uint32_t frameIndex_;
};

I3_MODULE(I3RNGKeyInjector);
