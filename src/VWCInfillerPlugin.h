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
        void prepareInfillForSlice(int id, const AAB<2, int> &xy_slice_box, float layer_height_mm, double layer_thickness_mm, int brush) override;

        // called to generate the infill in a surface
        bool generateInfill(int id, float layer_height_mm, double layer_thickness_mm, int brush, const ClipperLib::Paths &surface, std::vector<std::unique_ptr<IceSLInterface::IPath> > &fills,bool &preserve_order) override;

    protected:

        // various parameters
        bool verbose_ = false;

        // work parameters
        int extruder_id_ = -1;
        int num_beads_ = 0;
        double nozzle_diameter_ = 0.4; // mm
        double xy_mm_per_unit_ = 1.0f/4096.0f; // mm
        double minBeadWidth_ = 0.75 * nozzle_diameter_; // mm
        double maxBeadWidth_ = 2.0 * nozzle_diameter_; // mm
        const double simplificationThreshold_ = 1.05; // unitless;
};

class VWCInfillerPlugin : public IceSLInterface::IInfillerPlugin
{
    private:

        std::vector<float> minBeadWidth_;
        float minBeadWidth_min_ = 0.2f;
        float minBeadWidth_max_ = 0.4f;
        std::vector<float> maxBeadWidth_;
        float maxBeadWidth_min_ = 0.6f;
        float maxBeadWidth_max_ = 1.2f;
        std::vector<int> numBeads_;
        int minBeads_ = 0;
        int maxBeads_ = 5000;

    public:

        using IInfillerInterface = IceSLInterface::IInfillerInterface;

        VWCInfillerPlugin() = default;
        virtual ~VWCInfillerPlugin() = default;

        bool initialize(IceSLInterface::IPluginEnvironment &env) override {
            std::cerr << "VWC Infiller plugin initializing...\n";
            IInfillerPlugin::initialize(env);
            ImGui::SetCurrentContext(env.getImGuiContext());
            gluxInit();
            std::cerr << "VWC Infiller plugin initialized\n";
            return true;
        }
        void dispose() override {}
        void gui(bool postService) override {}

        std::string name()    const override { return "Variable-width contouring"; }
        std::string author()  const override { return "Samuel Hornus and Tim Kuipers"; }
        std::string comment() const override { return "See article \"Variable-width contouring for additive manufacturing\""; }
        std::string guid()    const override { return ""; }

        bool addPluginSettings(IceSLInterface::EnumerableSettingsInterface& enumerable_settings) override;

        std::unique_ptr<IInfillerInterface> createInfiller () override
        {
            return std::unique_ptr<VWCInfiller>(new VWCInfiller());
        }
};

// -----------------------------------------------
