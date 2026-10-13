import java.nio.file.*;
import org.objectweb.asm.*;

/** Version-pinned repair: use the node's existing per-row native cleaner. */
public class RegisterReactionVector {
    static final String NODE = "org/rdkit/knime/nodes/onecomponentreaction2/AbstractRDKitReactionNodeModel";
    public static void main(String[] args) throws Exception {
        ClassReader reader = new ClassReader(Files.readAllBytes(Path.of(args[0])));
        ClassWriter writer = new ClassWriter(reader, ClassWriter.COMPUTE_MAXS);
        final int[] sites = {0, 0};
        reader.accept(new ClassVisitor(Opcodes.ASM9, writer) {
            @Override public MethodVisitor visitMethod(int access, String name, String descriptor,
                    String signature, String[] exceptions) {
                MethodVisitor parent = super.visitMethod(access, name, descriptor, signature, exceptions);
                if (name.equals("<clinit>")) {
                    return new MethodVisitor(Opcodes.ASM9, parent) {
                        @Override public void visitCode() {
                            super.visitCode();
                            super.visitFieldInsn(Opcodes.GETSTATIC, "java/lang/System", "err", "Ljava/io/PrintStream;");
                            super.visitLdcInsn("GLADE_VECTOR_CLEANUP_PATCH_LOADED v1");
                            super.visitMethodInsn(Opcodes.INVOKEVIRTUAL, "java/io/PrintStream", "println", "(Ljava/lang/String;)V", false);
                            sites[1]++;
                        }
                    };
                }
                if (!name.equals("processReactionResults")) return parent;
                return new MethodVisitor(Opcodes.ASM9, parent) {
                    @Override public void visitMethodInsn(int opcode, String owner, String method,
                            String desc, boolean isInterface) {
                        super.visitMethodInsn(opcode, owner, method, desc, isInterface);
                        if (owner.equals("org/RDKit/ChemicalReaction") && method.equals("runReactants")
                                && desc.equals("(Lorg/RDKit/ROMol_Vect;)Lorg/RDKit/ROMol_Vect_Vect;")) {
                            super.visitVarInsn(Opcodes.ALOAD, 0);
                            super.visitInsn(Opcodes.SWAP);
                            super.visitVarInsn(Opcodes.LLOAD, 5);
                            super.visitMethodInsn(Opcodes.INVOKEVIRTUAL, NODE, "markForCleanup",
                                    "(Ljava/lang/Object;J)Ljava/lang/Object;", false);
                            super.visitTypeInsn(Opcodes.CHECKCAST, "org/RDKit/ROMol_Vect_Vect");
                            sites[0]++;
                        }
                    }
                };
            }
        }, 0);
        if (sites[0] != 1 || sites[1] != 1)
            throw new IllegalStateException("Unexpected bytecode: cleanup=" + sites[0] + ", init=" + sites[1]);
        Files.write(Path.of(args[1]), writer.toByteArray());
        System.out.println("PATCHED_CALL_SITES=" + sites[0] + " LOAD_MARKERS=" + sites[1]);
    }
}
