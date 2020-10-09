#pragma once

#include "IPlugin.h"
#include "IInfillerPlugin.h"

// -----------------------------------------------

extern "C" ICESL_PLUGIN_API IceSLInterface::IPlugin * createPlugin();

extern "C" ICESL_PLUGIN_API void destroyPlugin(IceSLInterface::IPlugin* plugin);

// -----------------------------------------------

class VWCInfiller : public IceSLInterface::IInfillerInterface
{
    public:

        VWCInfiller() = default;
        virtual ~VWCInfiller() = default;

        void setNumSlices(int num) override;
        void setWorldSpace(v3f world_corner_mm, v3f extent_mm) override;

        // called once at start
        void prepareInfill(int brush) override;

        // called once at end
        void terminateInfill(int brush) override;

        // called before each new slice
        void prepareInfillForSlice(int id, const AAB<2, int> &xy_slice_box, float height, int brush) override;

        // called to generate the infill in a surface
        bool generateInfill(int id, float slice_height_mm, int brush, const ClipperLib::Paths &surface, std::vector<std::unique_ptr<IceSLInterface::IPath> > &fills,bool &preserve_order) override;

    protected:

        // various parameters
        bool debug_ = false;
        bool verbose_ = false;

        // work parameters
        double nozzle_diameter_ = 0.4; // mm
        double xy_mm_per_unit_ = 1.0f/4096.0f; // mm
        double minBeadWidth_ = 0.75 * nozzle_diameter_; // mm
        double maxBeadWidth_ = 2.5 * nozzle_diameter_; // mm
        const double simplificationThreshold_ = 1.05; // unit-less;
};

class VWCInfillerPlugin : public IceSLInterface::IInfillerPlugin
{
    private:

        //float m_LineWidth_min_mm = 0.3f;
        //float m_LineWidth_max_mm = 1.0f;

    public:

        using IInfillerInterface = IceSLInterface::IInfillerInterface;

        VWCInfillerPlugin() = default;
        virtual ~VWCInfillerPlugin() = default;

        bool initialize(IceSLInterface::IPluginEnvironment &env) override {
            std::cerr << "VWC Infiller plugin initializing...\n";
            IInfillerPlugin::initialize(env);
            gluxInit();
            std::cerr << "VWC Infiller plugin initialized\n";
            return true;
        }
        void dispose() override {}

        std::string name()    const override { return "Variable-width Contouring"; }
        std::string author()  const override { return "Samuel Hornus and Tim Kuipers"; }
        std::string comment() const override { return "See article \"Variable-width contouring for additive manufacturing\""; }

        bool addInfillerSettings(IceSLInterface::EnumerableSettingsInterface& enumerable_settings) override;

        std::unique_ptr<IInfillerInterface> createInfiller () const override
        {
            return std::unique_ptr<VWCInfiller>(new VWCInfiller());
        }
};

// -----------------------------------------------
