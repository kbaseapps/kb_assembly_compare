
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
 * <p>Original spec-file type: Benchmark_assemblies_against_genomes_with_MUMmer4_Params</p>
 * <pre>
 * benchmark_assemblies_against_genomes_with_MUMmer4()
 * **
 * **  Align benchmark genomes to assembly contigs
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_genome_refs",
    "input_assembly_refs",
    "desc"
})
public class BenchmarkAssembliesAgainstGenomesWithMUMmer4Params {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_genome_refs")
    private String inputGenomeRefs;
    @JsonProperty("input_assembly_refs")
    private String inputAssemblyRefs;
    @JsonProperty("desc")
    private String desc;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public BenchmarkAssembliesAgainstGenomesWithMUMmer4Params withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_genome_refs")
    public String getInputGenomeRefs() {
        return inputGenomeRefs;
    }

    @JsonProperty("input_genome_refs")
    public void setInputGenomeRefs(String inputGenomeRefs) {
        this.inputGenomeRefs = inputGenomeRefs;
    }

    public BenchmarkAssembliesAgainstGenomesWithMUMmer4Params withInputGenomeRefs(String inputGenomeRefs) {
        this.inputGenomeRefs = inputGenomeRefs;
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

    public BenchmarkAssembliesAgainstGenomesWithMUMmer4Params withInputAssemblyRefs(String inputAssemblyRefs) {
        this.inputAssemblyRefs = inputAssemblyRefs;
        return this;
    }

    @JsonProperty("desc")
    public String getDesc() {
        return desc;
    }

    @JsonProperty("desc")
    public void setDesc(String desc) {
        this.desc = desc;
    }

    public BenchmarkAssembliesAgainstGenomesWithMUMmer4Params withDesc(String desc) {
        this.desc = desc;
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
        return ((((((((((("BenchmarkAssembliesAgainstGenomesWithMUMmer4Params"+" [workspaceName=")+ workspaceName)+", inputGenomeRefs=")+ inputGenomeRefs)+", inputAssemblyRefs=")+ inputAssemblyRefs)+", desc=")+ desc)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
