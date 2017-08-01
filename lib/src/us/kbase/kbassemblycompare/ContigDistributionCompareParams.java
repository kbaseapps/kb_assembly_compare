
package us.kbase.kbassemblycompare;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: Contig_Distribution_Compare_Params</p>
 * <pre>
 * contig_distribution_compare()
 * **
 * **  Compare Assembly Contig Length Distributions
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_assembly_refs"
})
public class ContigDistributionCompareParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_assembly_refs")
    private String inputAssemblyRefs;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public ContigDistributionCompareParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_assembly_refs")
    public String getInputAssemblyRefs() {
        return inputAssemblyRefs;
    }

    @JsonProperty("input_assembly_refs")
    public void setInputAssemblyRefs(String inputAssemblyRefs) {
        this.inputAssemblyRefs = inputAssemblyRefs;
    }

    public ContigDistributionCompareParams withInputAssemblyRefs(String inputAssemblyRefs) {
        this.inputAssemblyRefs = inputAssemblyRefs;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((("ContigDistributionCompareParams"+" [workspaceName=")+ workspaceName)+", inputAssemblyRefs=")+ inputAssemblyRefs)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
