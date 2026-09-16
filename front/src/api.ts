const errorMessage = (value: unknown) =>
  typeof value === "object" && value && "error" in value
    ? String(
        (value as { error: { message?: string } }).error.message || "请求失败",
      )
    : "请求失败";

export async function api<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(path, {
    headers: { "Content-Type": "application/json" },
    ...init,
  });
  const body = await response.json();
  if (!response.ok) {
    const error = new Error(errorMessage(body)) as Error & { status?: number };
    error.status = response.status;
    throw error;
  }
  return body as T;
}
