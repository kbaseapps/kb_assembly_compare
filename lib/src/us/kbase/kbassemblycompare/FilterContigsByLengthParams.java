
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
 * <p>Original spec-file type: Filter_Contigs_by_Length_Params</p>
 * <pre>
 * filter_contigs_by_length()
 * **
 * **  Remove Contigs that are under a minimum threshold
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_assembly_refs",
    "min_contig_length",
    "output_name"
})
public class FilterContigsByLengthParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_assembly_refs")
    private String inputAssemblyRefs;
    @JsonProperty("min_contig_length")
    private Long minContigLength;
    @JsonProperty("output_name")
    private String outputName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public FilterContigsByLengthParams withWorkspaceName(String workspaceName) {
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

    public FilterContigsByLengthParams withInputAssemblyRefs(String inputAssemblyRefs) {
        this.inputAssemblyRefs = inputAssemblyRefs;
        return this;
    }

    @JsonProperty("min_contig_length")
    public Long getMinContigLength() {
        return minContigLength;
    }

    @JsonProperty("min_contig_length")
    public void setMinContigLength(Long minContigLength) {
        this.minContigLength = minContigLength;
    }

    public FilterContigsByLengthParams withMinContigLength(Long minContigLength) {
        this.minContigLength = minContigLength;
        return this;
    }

    @JsonProperty("output_name")
    public String getOutputName() {
        return outputName;
    }

    @JsonProperty("output_name")
    public void setOutputName(String outputName) {
        this.outputName = outputName;
    }

    public FilterContigsByLengthParams withOutputName(String outputName) {
        this.outputName = outputName;
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
        return ((((((((((("FilterContigsByLengthParams"+" [workspaceName=")+ workspaceName)+", inputAssemblyRefs=")+ inputAssemblyRefs)+", minContigLength=")+ minContigLength)+", outputName=")+ outputName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
